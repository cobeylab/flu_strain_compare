from pymol import cmd, movie
from os.path import exists
from Bio import SeqIO
import pandas as pd
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
        is_aligned,
        position_map):
        assert (exists(query_sequence_file)), "Query sequence file must exist"
        assert (type(lineage) == str), "Lineage must be a string."
        assert (type(name) == str), "Name must be a string."
        assert (type(query_sequence_id) == str), "Query sequence ID must be a string."
        assert (type(is_aligned) == bool), "is_aligned must be a boolean."
        query_seqs = {s.id: s for s in SeqIO.parse(query_sequence_file, "fasta")}
        assert (len(query_seqs) >= 1), "Query sequence file must contain at least 1 sequence."
        self.name = name
        self.lineage = lineage
        self.sequence = query_seqs[query_sequence_id]
        self.is_aligned = is_aligned
        self.position_map = position_map
        self.query_sequence_file = query_sequence_file

        # Get PNGS sites
        gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
        self.pngs = [position_map.loc[m.start(), "H3"] for m in gly.finditer(str(self.sequence.seq))]

class seq_compare:
    def __init__(self, seq1, seq2):
        assert (type(seq1) == flu_seq), "Seq1 must be a flu_seq object"
        assert (type(seq2) == flu_seq), "Seq2 must be a flu_seq object"
        assert (seq1.is_aligned), "Seq1 must be aligned"
        assert (seq2.is_aligned), "Seq2 must be aligned"
        assert (seq1.query_sequence_file == seq2.query_sequence_file), "Sequences must be from the same alignment"
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

position_map_infile = "../data/EPI1487157_H3_Conversion.txt"
aligned_seqs = "../data/H3N2_HA_GISAID.fasta"
q1_id = "EPI1487157"
q1_name = "A/Minnesota/41/2019"
q2_id = "EPI1752480"
q2_name = "A/Tasmania/503/2020"
seq_lineage = "H3N2"

position_map = pd.read_csv(position_map_infile, sep = "\t")
position_map = position_map[position_map.Query != "-"].reset_index()

s1 = flu_seq(name = q1_name,
    lineage = seq_lineage,
    query_sequence_file = aligned_seqs,
    query_sequence_id = q1_id,
    is_aligned = True,
    position_map = position_map
    )

s2 = flu_seq(name = q2_name,
    lineage = seq_lineage,
    query_sequence_file = aligned_seqs,
    query_sequence_id = q2_id,
    is_aligned = True,
    position_map = position_map
    )

mutation_list = seq_compare(seq1 = s1, seq2 = s2).identify_mutations()
mutations = [m.position for m in mutation_list if m.position != "-"]
glycosylations = seq_compare(seq1 = s1, seq2 = s2).identify_PNGS_changes()

cmd.reinitialize()
cmd.fetch('4we8')
cmd.symexp('sym','4we8','4we8','3')
cmd.delete('sym11000000')
cmd.delete('sym04000000')
cmd.center('4we8, sym01000000, sym02000000')
cmd.set('ray_trace_mode', 0)
cmd.hide('everything')
cmd.show('surface')
cmd.remove('solvent')
cmd.color('gray70', 'all')
cmd.set('label_shadow_mode', 2)
cmd.set('label_position', (0, 0, 20))
cmd.set('label_size', -2)
cmd.set('label_color', 'black')

cmd.select('mutations', '(resi %s)'%'+'.join([i for i in mutations]))
cmd.select('glycan_additions', '(resi %s)'%'+'.join([i for i in glycosylations['glycan_additions']]))
cmd.select('glycans_shared', '(resi %s)'%'+'.join([i for i in glycosylations['glycans_shared']]))


cmd.color('blue', 'glycans_shared')
cmd.color('green', 'glycan_additions')

cmd.color('yellow', 'mutations')

for m in mutation_list:
    label = m.mutation
    resi = m.position
    if resi != "-":
        cmd.select(label, 'n. CA and i. ' + resi)
        label_name = str(label)
        cmd.label(selection = label, expression = "label_name")

for g in glycosylations['glycan_additions']:
    label = "PNGS" + str(g)
    resi = str(g)
    if (resi != "-"):
        cmd.select(label, 'n. CA and i. ' + resi)
        label_name = str(label)
        cmd.label(selection = label, expression = "label_name")

for g in glycosylations['glycans_shared']:
    label = "PNGS" + str(g)
    resi = str(g)
    if (resi != "-"):
        cmd.select(label, 'n. CA and i. ' + resi)
        label_name = str(label)
        cmd.label(selection = label, expression = "label_name")
cmd.set_view((
    -0.424973905,   -0.644264817,    0.635839999,
    -0.903248310,    0.255860180,   -0.344455302,
     0.059236884,   -0.720711648,   -0.690678000,
     0.000000000,    0.000000000, -417.670715332,
     0.000007629,  -58.260734558,   10.906387329,
   329.294769287,  506.046661377,  -20.000000000 ))
cmd.png('../figures/%s-%s.png'%(q1_name.replace("/","_"), q2_name.replace("/","_")), ray=1, dpi=300)