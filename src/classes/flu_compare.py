from os.path import exists
from Bio import SeqIO
import re
from os import system
import pandas as pd

class flu_mutation:
    def __init__(self,
        pymol_resi,
        label,
        strain1,
        strain2):
        assert (type(pymol_resi) == str), "Position must be a string"
        assert (type(label) == str), "Mutation must be identified as a string"
        assert (type(strain1) == str), "Strain 1 must be a string identifier"
        assert (type(strain2) == str), "Strain 2 must be a string identifier"
        self.pymol_resi = pymol_resi
        self.label = label
        self.strain1 = strain1
        self.strain2 = strain2
    def __str__(self):
        return f"Strain 1: {self.strain1}, Strain 2: {self.strain2}, Mutation: {self.label}"

class flu_pngs:
    def __init__(self,
        pymol_resi,
        label):
        assert (type(pymol_resi) == str), "Position must be a string"
        assert (type(label) == str), "Label must be identified as a string"
        self.pymol_resi = pymol_resi
        self.label = label


class flu_seq:
    def __init__(self,
        name,
        lineage,
        query_sequence_file,
        query_sequence_id):
        assert (exists(query_sequence_file)), "Query sequence file must exist"
        assert (type(lineage) == str), "Lineage must be a string."
        assert (type(name) == str), "Name must be a string."
        assert (type(query_sequence_id) == str), "Query sequence ID must be a string."
        query_seqs = {s.id: s for s in SeqIO.parse(query_sequence_file, "fasta")}
        assert (len(query_seqs) >= 1), "Query sequence file must contain at least 1 sequence."
        self.name = name
        self.lineage = lineage
        self.sequence = query_seqs[query_sequence_id]
        self.query_sequence_file = query_sequence_file
        self.align_to_reference()

        # Get PNGS sites
        gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
        self.pngs = [str(m.start() + 1) + "_" for m in gly.finditer(str(self.sequence.seq))]
        
    def align_to_reference(self):
        ref_file = "/app/data/%s_ref.fasta"%(self.lineage)
        temp_seqfile = "/app/figures/tmp.fasta"
        temp_alignfile = "/app/figures/aligned.fasta"
        SeqIO.write([self.sequence], temp_seqfile, "fasta")
        command = "mafft --keeplength --add %s %s > %s"%(temp_seqfile, ref_file, temp_alignfile)
        system(command)
        newseq = [s for s in SeqIO.parse(temp_alignfile, "fasta") if s.id == self.sequence.id]
        assert (len(newseq) == 1), "Alignment not found"
        self.sequence = newseq[0]
        
class seq_compare:
    def __init__(self, seq1, seq2, numbering_scheme):
        assert (type(seq1) == flu_seq), "Seq1 must be a flu_seq object"
        assert (type(seq2) == flu_seq), "Seq2 must be a flu_seq object"
        assert (seq1.lineage == seq2.lineage), "Cannot compare sequences from different lineages"
        self.seq1 = seq1
        self.seq2 = seq2
        self.lineage = seq1.lineage
        self.numbering_scheme = numbering_scheme
        assert (exists("/app/data/%s_Conversion.csv"%self.lineage)), "Conversion file does not exist"
        self.conversion_table = pd.read_csv("/app/data/%s_Conversion.csv"%self.lineage,
            index_col = "ref_one_index")
    def convert_numbering(self,
        position):
        return self.conversion_table.loc[position + 1, self.numbering_scheme]

    def identify_mutations(self):
        mutations_out = []
        for i, (b1, b2) in enumerate(zip(self.seq1.sequence.seq, self.seq2.sequence.seq)):
            if b1 != b2:
                p = self.convert_numbering(i)
                mutations_out.append(
                    flu_mutation(pymol_resi = str(i+1) + "_",
                        label = "".join([str(b1), str(p), str(b2)]),
                        strain1 = self.seq1.name,
                        strain2 = self.seq2.name)
                )
        return mutations_out
    def identify_PNGS_changes(self, change_type):
        if change_type == "deletions":
            comparison_set = set(self.seq1.pngs) - set(self.seq2.pngs)
        elif change_type == "additions":
            comparison_set = set(self.seq2.pngs) - set(self.seq1.pngs)
        elif change_type == "shared":
            comparison_set = set(self.seq2.pngs).intersection(set(self.seq1.pngs))
        pngs_out = []
        for pymol_position in comparison_set:
            conversion_index = int(pymol_position.replace("_","")) - 1
            label = "PNGS%s"%self.convert_numbering(conversion_index)
            pngs_out.append(flu_pngs(pymol_resi = pymol_position,
                label = label))
        return pngs_out
