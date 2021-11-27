from os.path import exists
from Bio import SeqIO
import re
from os import system
import pandas as pd
from pymol import cmd

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
    
    def make_figure(self):
        q1_name = self.seq1.name
        q2_name = self.seq2.name
        mutation_list = self.identify_mutations()
        mutations = [m.pymol_resi for m in mutation_list if m.pymol_resi != "-"]
        gly_del = self.identify_PNGS_changes("deletions")
        gly_add = self.identify_PNGS_changes("additions")
        gly_share = self.identify_PNGS_changes("shared")

        cmd.load('/app/data/%s_renumbered.pse'%self.lineage)
        cmd.set('ray_trace_mode', 0)

        # Label parameters
        cmd.set('label_position', (0, 0, 20))
        cmd.set('label_size', -4)
        cmd.set('label_color', 'black')

        # Color mutations
        cmd.select('mutations', '(resi %s)'%'+'.join([i for i in mutations]))
        cmd.color('yellow', 'mutations')

        # Color glycosylations
        color_pngs(gly_del, "glycan_deletions", "red")
        color_pngs(gly_add, "glycan_additions", "green")
        color_pngs(gly_share, "glycans_shared", "blue")

        # Add labels
        label_resi(mutation_list)
        label_resi(gly_del)
        label_resi(gly_add)
        label_resi(gly_share)

        y_start = -60
        x_pos = -50
        z_pos = 10
        create_label(x = x_pos,
            y = y_start,
            z = z_pos,
            label_text = "Mutations",
            label_name = "mutation_label",
            label_color = "yellow")
        create_label(x = x_pos,
            y = y_start - 3,
            z = z_pos,
            label_text = "Shared PNGS",
            label_name = "pngs_shared_label",
            label_color = "blue")
        create_label(x = x_pos,
            y = y_start - 6,
            z = z_pos,
            label_text = "Added PNGS",
            label_name = "pngs_add_label",
            label_color = "green")
        create_label(x = x_pos,
            y = y_start - 9,
            z = z_pos,
            label_text = "Deleted PNGS",
            label_name = "pngs_del_label",
            label_color = "red")

        base_filename = "%s-%s"%(q1_name.replace("/","_"),
                q2_name.replace("/","_"))
        return base_filename

def color_pngs(glylist, name, color):
    if len(glylist) > 0:
        cmd.select(name, '(resi %s)'%'+'.join([g.pymol_resi for g in glylist]))
        cmd.color(color, name)
def create_label(x, y, z, label_text, label_name, label_color):
    cmd.pseudoatom(label_name, pos=[x, y, z])
    global label_temp
    label_temp = label_text
    cmd.label(selection = label_name, expression = "label_temp")
    cmd.hide("wire", selection = label_name)
    cmd.set("label_color", selection = label_name, value = label_color)
def label_resi(resilist):
    for m in resilist:
        label = m.label
        resi = m.pymol_resi
        if resi != "-":
            cmd.select(label, 'n. CA and i. ' + resi)
            global label_name
            label_name = str(label)
            cmd.label(selection = label, expression = "label_name")
def make_comparison_object(parameters):
    seq_file = "/app/data/" + parameters["seq_file"]
    q1_id = parameters["q1_id"]
    q1_name = parameters["q1_name"]
    q2_id = parameters["q2_id"]
    q2_name = parameters["q2_name"]
    seq_lineage = parameters["seq_lineage"]
    numbering_scheme = parameters["numbering_scheme"]
    s1 = flu_seq(name = q1_name,
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q1_id
        )
    s2 = flu_seq(name = q2_name,
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q2_id
        )
    comparison = seq_compare(seq1 = s1,
        seq2 = s2,
        numbering_scheme = numbering_scheme)
    return comparison
