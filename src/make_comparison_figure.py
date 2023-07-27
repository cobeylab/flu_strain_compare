from pymol.viewing import color
from pymol import cmd
from classes.flu_compare import make_comparison_object, make_figure, SequenceComparisonEncoder
import json
import sys

parameters = json.load(open(sys.argv[1]))
figure_dir = parameters["output_dir"]

comparison = make_comparison_object(parameters)

# Output json or csv
if parameters.get("output_csv", False):
    print("PyMOL residue\tMutation\tPercent conserved\tConservative")
    [print(m.as_row()) for m in comparison.mutation_list]
else:
    print(json.dumps(comparison, indent=4, cls=SequenceComparisonEncoder))


base_filename = make_figure(comparison)

# Save figures
cmd.set("opaque_background", "on")

exports = parameters["export_files"]

if "PNG" in (e.upper() for e in exports):
    cmd.png('%s/%s.png'%(
            figure_dir,
            base_filename),
        width=800,
        height=800,
        ray=1,
        dpi=300)

if "PSE" in (e.upper() for e in exports):
    cmd.save('%s/%s.pse'%(
            figure_dir,
            base_filename))
