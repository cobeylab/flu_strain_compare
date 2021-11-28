import sys
sys.path.append("/app/src/classes")
from os import system
from pymol.viewing import color
from pymol import cmd
from flu_compare import make_comparison_object, flu_seq
import label_globals
import json

parameters = json.load(open("/app/configuration/movie_config.json"))
all_pngs = []

system("mkdir -p frames")

for i, s_id in enumerate(parameters["frame_order"]):
    if i != len(parameters["frame_order"]) - 1:
        parameters["q1_id"] =  s_id
        parameters["q2_id"] = parameters["frame_order"][i + 1]
        comparison = make_comparison_object(parameters)
        base_filename = comparison.make_figure()
        # Save figures
        cmd.png('frames/%s_%s.png'%(
                str(i+1),
                base_filename),
            width=800,
            height=800,
            ray=1,
            dpi=300)

system("convert -dispose none -delay 0 -size 800x800 xc:Black +antialias -delay 100 -loop 0 -dispose previous -layers OptimizePlus frames/*.png /app/figures/%s.mp4"%parameters["output_handle"])
system("rm frames/*.png")