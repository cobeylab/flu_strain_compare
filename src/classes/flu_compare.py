from os.path import exists
from Bio import SeqIO
import re

class flu_mutation:
    def __init__(self,
        position,
        mutation,
        strain1,
        strain2):
        assert (type(position) == str), "Position must be a string"
        assert (type(mutation) == str), "Mutation must be identified as a string"
        assert (type(strain1) == str), "Strain 1 must be a string identifier"
        assert (type(strain2) == str), "Strain 2 must be a string identifier"
        self.position = position
        self.mutation = mutation
        self.strain1 = strain1
        self.strain2 = strain2
    def __str__(self):
        return f"Strain 1: {self.strain1}, Strain 2: {self.strain2}, Mutation: {self.mutation}"

class flu_seq:
    def __init__(self,
        name,
        lineage,
        query_sequence_file,
        query_sequence_id,
        position_map):
        assert (exists(query_sequence_file)), "Query sequence file must exist"
        assert (type(lineage) == str), "Lineage must be a string."
        assert (type(name) == str), "Name must be a string."
        assert (type(query_sequence_id) == str), "Query sequence ID must be a string."
        query_seqs = {s.id: s for s in SeqIO.parse(query_sequence_file, "fasta")}
        assert (len(query_seqs) >= 1), "Query sequence file must contain at least 1 sequence."
        assert ("H3" in position_map.columns and "Query" in position_map.columns), "Position map must contain columns named H3 and Query"
        self.name = name
        self.lineage = lineage
        self.sequence = query_seqs[query_sequence_id]
        assert (len(self.sequence) == 566), "For now, sequence length must be 566."
        self.position_map = position_map
        self.query_sequence_file = query_sequence_file

        # Get PNGS sites
        gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
        self.pngs = [position_map.loc[m.start(), "H3"] for m in gly.finditer(str(self.sequence.seq))]
        
    def align_to_reference(self):
        ref_file = "%s_ref.fasta"%(self.lineage)
        temp_seqfile = "tmp/tmp.fasta"
        temp_alignfile = "tmp/aligned.fasta"
        command = "mafft --keeplength %s %s > %s"(ref_file, temp_seqfile, temp_alignfile)
        newseq = [s for s in SeqIO.parse(temp_alignfile, "fasta") if s.id == self.sequence.id]
        assert (len(newseq) == 1), "Alignment not found"
        self.sequence = newseq[0]
        
class seq_compare:
    def __init__(self, seq1, seq2):
        assert (type(seq1) == flu_seq), "Seq1 must be a flu_seq object"
        assert (type(seq2) == flu_seq), "Seq2 must be a flu_seq object"
        self.seq1 = seq1
        self.seq2 = seq2
    def identify_mutations(self):
        mutations_out = []
        for i, (b1, b2) in enumerate(zip(self.seq1.sequence.seq, self.seq2.sequence.seq)):
            if b1 != b2:
                assert (self.seq1.position_map.loc[i, "H3"] == self.seq2.position_map.loc[i, "H3"]), "Sequences must use same position map"
                p = self.seq1.position_map.loc[i, "H3"]
                mutations_out.append(
                    flu_mutation(position = str(p),
                        mutation = "".join([str(b1), str(p), str(b2)]),
                        strain1 = self.seq1.name,
                        strain2 = self.seq2.name)
                )
        return mutations_out
    def identify_PNGS_changes(self):
        glycan_deletions = set(self.seq1.pngs) - set(self.seq2.pngs)
        glycan_additions = set(self.seq2.pngs) - set(self.seq1.pngs)
        glycans_shared = set(self.seq2.pngs).intersection(set(self.seq1.pngs))
        return {"glycan_deletions": glycan_deletions, "glycan_additions": glycan_additions, "glycans_shared": glycans_shared}  
