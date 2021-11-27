import sys
sys.path.append("/app/src/classes")
from pymol.viewing import color
from pymol import cmd
from utilities import *
import json

parameters = json.load(open("/app/configuration/config.json"))
figure_dir = "/app/figures/"

mutations, gly_del, gly_add, gly_share, mutation_list, q1_name, q2_name, seq_lineage = make_comparison_object(parameters)

cmd.load('/app/data/%s_renumbered.pse'%seq_lineage)
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
    label_text="Deleted PNGS",
    label_name="pngs_del_label",
    label_color="red")

# Save figures
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