
import sys
sys.path.append("/app/src/classes")
from pymol.viewing import color
from pymol import cmd
from flu_compare import make_comparison_object, label_resi, create_label
import label_globals
import json

parameters = json.load(open("/app/configuration/config.json"))
figure_dir = "/app/figures/"

comparison = make_comparison_object(parameters)
base_filename = comparison.make_figure()

# Save figures
cmd.png('%s/%s.png'%(
        figure_dir,
        base_filename),
    width=800,
    height=800,
    ray=1,
    dpi=300)
cmd.save('%s/%s.pse'%(
        figure_dir,
        base_filename))