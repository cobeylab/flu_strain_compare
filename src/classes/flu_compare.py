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
        lineage,
        query_sequence_file,
        query_sequence_id):
        assert (exists(query_sequence_file)), "Query sequence file must exist"
        assert (type(lineage) == str), "Lineage must be a string."
        assert (type(query_sequence_id) == str), "Query sequence ID must be a string."
        query_seqs = {s.id: s for s in SeqIO.parse(query_sequence_file, "fasta")}
        assert (len(query_seqs) >= 1), "Query sequence file must contain at least 1 sequence."
        self.lineage = lineage
        self.sequence = query_seqs[query_sequence_id]
        self.name = self.sequence.description.split(" | ")[1]
        self.query_sequence_file = query_sequence_file
        self.align_to_reference()

        # Get PNGS sites
        gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
        self.pngs = [str(m.start() + 1) for m in gly.finditer(str(self.sequence.seq))]
        
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
        self.mutation_list = self.identify_mutations()
        self.gly_del = self.identify_PNGS_changes("deletions")
        self.gly_add = self.identify_PNGS_changes("additions")
        self.gly_share = self.identify_PNGS_changes("shared")
    def convert_numbering(self,
        position):
        return self.conversion_table.loc[position + 1, self.numbering_scheme]

    def identify_mutations(self):
        mutations_out = []
        for i, (b1, b2) in enumerate(zip(self.seq1.sequence.seq, self.seq2.sequence.seq)):
            if b1 != b2:
                p = self.convert_numbering(i)
                mutations_out.append(
                    flu_mutation(pymol_resi = str(i+1),
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
        
        mutations = [m.pymol_resi for m in self.mutation_list if m.pymol_resi != "-"]
        cmd.reinitialize()
        cmd.set('ray_trace_mode', 0)
        # Label parameters
        cmd.set('label_position', (0, 0, 20))
        cmd.set('label_size', -3)
        cmd.set('label_color', 'black')
        cmd.load('/app/data/%s_pngs.pse'%self.lineage)

        # Color mutations
        cmd.select('mutations', '(resi %s)'%'+'.join([i for i in mutations]))
        cmd.color('yellow', 'mutations')

        # Color glycosylations
        color_pngs(self.gly_del, "glycan_deletions", "red")
        color_pngs(self.gly_add, "glycan_additions", "green")
        color_pngs(self.gly_share, "glycans_shared", "blue")
        base_filename = "%s-%s"%(q1_name.replace("/","_"),
                q2_name.replace("/","_"))
        # Add labels
        label_resi(self.mutation_list)
        label_resi(self.gly_del)
        label_resi(self.gly_add)
        label_resi(self.gly_share)

        x_pos = -60
        z_pos = 0
        y_start = 0
        y_offset = -5
        create_label(x_pos, y_start, z_pos, "Mutations", "mutlabel", "yellow")
        create_label(x_pos, y_start + y_offset, z_pos, "PNGS added", "glyaddlabel", "green")
        create_label(x_pos, y_start + 2*y_offset , z_pos, "PNGS deleted", "glydellabel", "red")
        create_label(x_pos, y_start + 3*y_offset, z_pos, "PNGS shared", "glysharelabel", "blue")
        cmd.hide("everything", "extra_glycans")

        create_label(0, 110, 0, "%s vs. %s"%(q1_name, q2_name), "strains", "white", label_size=-5)
        return base_filename

def color_pngs(glylist, name, color):
    if len(glylist) > 0:
        PNGS_names = ["PNGS%s"%g.pymol_resi for g in glylist]
        PNGS_names_final = set(PNGS_names).intersection(set(cmd.get_names(type="selections")))
        # Need to add warning here if it finds a PNGS site that isn't in the structure
        cmd.select(name, ' | '.join(PNGS_names_final))
        cmd.show("sticks", name)
        cmd.color(color, name)
def make_comparison_object(parameters):
    seq_file = "/app/data/" + parameters["seq_file"]
    q1_id = parameters["q1_id"]
    q2_id = parameters["q2_id"]
    seq_lineage = parameters["seq_lineage"]
    numbering_scheme = parameters["numbering_scheme"]
    s1 = flu_seq(
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q1_id
        )
    s2 = flu_seq(
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q2_id
        )
    comparison = seq_compare(seq1 = s1,
        seq2 = s2,
        numbering_scheme = numbering_scheme)
    return comparison

def create_label(x, y, z, label_text, label_name, label_color, label_size=-4):
    cmd.pseudoatom(label_name)
    cmd.label(selection = label_name, expression = f"'{label_text}'")
    cmd.hide("wire", selection = label_name)
    cmd.set("label_color", selection = label_name, value = label_color)
    cmd.set("label_position", selection = label_name, value = [x,y,z])
    cmd.set("label_size", selection = label_name, value = label_size)
def label_resi(resilist):
    for m in resilist:
        label = m.label
        resi = m.pymol_resi
        if resi != "-":
            cmd.select(label, 'n. CA and i. ' + resi)
            cmd.label(selection = label, expression = f"'{label}'")
