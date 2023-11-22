from pymol import cmd
import os
import sys

input_file_base = ""

# Parse parameters, open input file, derive output directory name.
try:
    input_param = sys.argv[1]
    print(f"Opening {input_param}...")
    cmd.load(input_param)
    input_file_base = os.path.split(input_param)[-1].split(".")[0]
except:
    print(f"Provide input .pse file path: python <path>/<filename>.pse.")
    exit()

# Create output directory based on input file name.
os.mkdir(input_file_base)

# Iterate over objects in file and print the type and properties of the atoms to files named for the PyMOL chain. Note: Does not include coordinates and a few other properties. Objects are generally instances of chempy classes.
# Chains are pymol models: https://fossies.org/dox/pymol-open-source-2.5.0/models_8py_source.html
for x in cmd.get_names("objects"):
    for ch in cmd.get_chains(x):
        m = cmd.get_model(f"chain {ch}")
        print(f"Saving objects in chain {ch}...")
        filename = os.path.join(input_file_base, f"chain_{ch}.txt")
        with open(filename, "a") as f:
            print(type(m), file=f)
            for at in m.atom:
                print(f"{at.resn} {at.resi} {at.name}", file=f)

# Iterate over selections and print to "selections.txt" file.
print("Saving selections...")
filename = os.path.join(input_file_base, f"selections.txt")
with open(filename, "a") as f:
    for x in cmd.get_names("selections"):
        m = cmd.get_model(x)
        print(x, file=f)
        print(type(m), file=f)
        for a in m.atom:
            print(f"{at.resn} {at.resi} {at.name}", file=f)

