from pymol.viewing import color
from pymol import cmd
from classes.flu_compare import make_comparison_object, SequenceComparisonEncoder
import json

parameters = json.load(open("configuration/config.json"))
figure_dir = "figures/"

comparison = make_comparison_object(parameters)

print(json.dumps(comparison, indent=4, cls=SequenceComparisonEncoder))

base_filename = comparison.make_figure()

# Save figures
cmd.set("opaque_background", "on")
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
