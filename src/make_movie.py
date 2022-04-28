import os
from pymol.viewing import color
from pymol import cmd
from classes.flu_compare import make_comparison_object, make_figure, SequenceComparisonEncoder
import json

parameters = json.load(open("configuration/movie_config.json"))
all_pngs = []

OUTPUT_DIR = parameters["output_dir"]
DATA_DIR = "data"
FRAMES_DIR = f"{DATA_DIR}/tmp/frames"

if not os.path.exists(FRAMES_DIR):
    os.makedirs(FRAMES_DIR)
pngs = []

for i, s_id in enumerate(parameters["frame_order"]):
    if i != len(parameters["frame_order"]) - 1:
        parameters["q1_id"] =  s_id
        parameters["q2_id"] = parameters["frame_order"][i + 1]
        comparison = make_comparison_object(parameters)
        base_filename = make_figure(comparison)
        png = f"{FRAMES_DIR}/{str(i+1)}_{base_filename}.png"
        pngs.append(png)
        # Save figures
        cmd.png(png,
            width=800,
            height=800,
            ray=1,
            dpi=300)

OUPUT_HANDLE = parameters["output_handle"]
convert_command = f"convert -dispose none -delay 0 -size 800x800 xc:Black +antialias -delay 100 -loop 0 -dispose previous -layers OptimizePlus {FRAMES_DIR}/*.png {OUTPUT_DIR}/{OUPUT_HANDLE}.mp4"

print(f"Converting frames:\n{convert_command}")
os.system(convert_command)

for frame_file in pngs:
    try:
        os.remove(frame_file)
    except:
        print(f"Warning: Unable to remove frame file {frame_file}")
os.rmdir(FRAMES_DIR)

