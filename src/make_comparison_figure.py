import sys
sys.path.append("./classes")
from pymol import cmd
from flu_compare import flu_seq,seq_compare
import pandas as pd

# Path to a csv file that maps positions in query 
# sequence to a numbering scheme, i.e., H3 numbering.
# You can easily generate such a file from the Influenza Research Database's 
# HA subtype numbering conversion tool.

position_map_infile = "../data/H3_Conversion.txt" 
seqs = "../data/H3_cell_vaccines.fasta"
q1_id = "EPI1409001"
q1_name = "A/Hong Kong/45/2019"
q2_id = "EPI1548699"
q2_name = "A/Minnesota/41/2019"
seq_lineage = "H3N2"

position_map = pd.read_csv(position_map_infile, sep = "\t")
position_map = position_map[position_map.Query != "-"].reset_index()

s1 = flu_seq(name = q1_name,
    lineage = seq_lineage,
    query_sequence_file = seqs,
    query_sequence_id = q1_id,
    position_map = position_map
    )

s2 = flu_seq(name = q2_name,
    lineage = seq_lineage,
    query_sequence_file = seqs,
    query_sequence_id = q2_id,
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
cmd.color('yellow', 'mutations')

if len(glycosylations['glycan_additions']) > 0:
    cmd.select('glycan_additions', '(resi %s)'%'+'.join([i for i in glycosylations['glycan_additions']]))
    cmd.color('green', 'glycan_additions')
if len(glycosylations['glycans_shared']) > 0:
    cmd.select('glycans_shared', '(resi %s)'%'+'.join([i for i in glycosylations['glycans_shared']]))
    cmd.color('blue', 'glycans_shared')

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