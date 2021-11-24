import sys

from pymol.viewing import color
sys.path.append("/usr/src/classes")
from pymol import cmd
from flu_compare import flu_seq,seq_compare
import pandas as pd
import json



parameters = json.load(open("/usr/configuration/config.json"))
figure_dir = "/usr/figures/"
seq_file = "/usr/data/" + parameters["seq_file"]
q1_id = parameters["q1_id"]
q1_name = parameters["q1_name"]
q2_id = parameters["q2_id"]
q2_name = parameters["q2_name"]
seq_lineage = parameters["seq_lineage"]

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

mutation_list = seq_compare(seq1 = s1, seq2 = s2).identify_mutations()
mutations = [m.pymol_resi for m in mutation_list if m.pymol_resi != "-"]
glycosylations = seq_compare(seq1 = s1, seq2 = s2).identify_PNGS_changes()

cmd.load('/usr/data/%s_renumbered.pse'%seq_lineage)
cmd.set('ray_trace_mode', 0)

# Label parameters
#cmd.set('label_shadow_mode', 2)
cmd.set('label_position', (0, 0, 20))
cmd.set('label_size', -4)
cmd.set('label_color', 'black')

# Color mutations
cmd.select('mutations', '(resi %s)'%'+'.join([i for i in mutations]))
cmd.color('yellow', 'mutations')

# Color glycosylations
if len(glycosylations['glycan_additions']) > 0:
    cmd.select('glycan_additions', '(resi %s)'%'+'.join([i for i in glycosylations['glycan_additions']]))
    cmd.color('green', 'glycan_additions')
if len(glycosylations['glycan_deletions']) > 0:
    cmd.select('glycan_deletions', '(resi %s)'%'+'.join([i for i in glycosylations['glycan_deletions']]))
    cmd.color('red', 'glycan_deletions')
if len(glycosylations['glycans_shared']) > 0:
    cmd.select('glycans_shared', '(resi %s)'%'+'.join([i for i in glycosylations['glycans_shared']]))
    cmd.color('blue', 'glycans_shared')

# Add in labels...should probably make some utility functions for this part
for m in mutation_list:
    label = m.mutation
    resi = m.pymol_resi
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
for g in glycosylations['glycan_deletions']:
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

def create_label(x, y, z, label_text, label_name, label_color):
    cmd.pseudoatom(label_name, pos=[x,y,z])
    global label_temp
    label_temp = label_text
    cmd.label(selection = label_name, expression = "label_temp")
    cmd.hide("wire", selection = label_name)
    cmd.set("label_color", selection = label_name, value = label_color)

y_start = -60

create_label(x=-50,
    y=y_start,
    z=10,
    label_text="Mutations",
    label_name="mutation_label",
    label_color="yellow")
create_label(x=-50,
    y=y_start - 3,
    z=10,
    label_text="Shared PNGS",
    label_name="pngs_shared_label",
    label_color="blue")
create_label(x=-50,
    y=y_start - 6,
    z=10,
    label_text="Added PNGS",
    label_name="pngs_add_label",
    label_color="green")
create_label(x=-50,
    y=y_start - 9,
    z=10,
    label_text="Deleted PNGS",
    label_name="pngs_del_label",
    label_color="red")


cmd.png('%s/%s-%s.png'%(
        figure_dir,
        q1_name.replace("/","_"),
        q2_name.replace("/","_")),
    width=800,
    height=800,
    ray=1,
    dpi=300)

cmd.save('%s/%s-%s.pse'%(
        figure_dir,
        q1_name.replace("/","_"),
        q2_name.replace("/","_")))