from pymol import cmd
import util
import re

cmd.reinitialize()
cmd.fetch('4M4Y')
cmd.remove('solvent')
cmd.hide('everything')
cmd.show('surface')
cmd.color('gray70', 'all')
cmd.rotate("y", "90")
cmd.rotate("z", "-45")
cmd.rotate("y", "15")
cmd.rotate("x", "15")

#cmd.remove("chain G+H")
#cmd.remove("resi 400+401+402+403+407")

myspace = {'resis': [], 'excluded_resi': ['400', '401', '402', '403', '407']}
resis = myspace['resis']
cmd.iterate("chain A", "if (resi, resn) not in resis and resi not in excluded_resi:\n\tresis.append((resi, resn))", space = myspace)
start_position = 14
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain A+C+E & resi %s & resn %s"%(r,aa), "resi='%s'"%(new_resi))


myspace = {'resis': [], 'excluded_resi': ['400', '401', '402', '403', '407']}
resis = myspace['resis']
cmd.iterate("chain B", "if (resi, resn) not in resis and resi not in excluded_resi:\n\tresis.append((resi, resn))", space = myspace)
start_position = 345
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain B+D+F & resi %s & resn %s"%(r,aa), "resi='%s'"%(str(i+start_position)))
cmd.iterate("chain B", "print(resi, resn)")


# Attempt to find and label PNGSes by regular expression.
seq_model = ""
m = cmd.get_model("chain E")
prev_index = None

# Assemble sequence.
for at in m.atom:
    if prev_index is None or at.resi != prev_index:
        if at.resn in util.AA3to1:
            seq_model += util.AA3to1[at.resn]
    prev_index = at.resi

print(seq_model)

# Find pattern and label selections.
gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
model_pngs = [str(m.start() + 1) for m in gly.finditer(str(seq_model))]

for p in model_pngs:
    try:
        idx = int(p) + 14
        cmd.select(f"PNGS{idx}", f"i. {idx}")
    except Exception as e:
        print(e)



# TODO kludge to avoid extra_glycans error
cmd.select("extra_glycans", "i. 0")

cmd.save("H1_minimal.pse")
