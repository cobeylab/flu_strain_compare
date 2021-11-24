from pymol import cmd
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

# Removing glycans
cmd.remove("chain G+H")
cmd.remove("resi 400+401+402+403+407")

resis = []
cmd.iterate("chain A", "if (resi, resn) not in resis:\n\tresis.append((resi, resn))")
start_position = 14
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain A+C+E & resi %s & resn %s"%(r,aa), "resi='%s_'"%(new_resi))
cmd.iterate("chain A", "print(resi, resn)")


resis = []
cmd.iterate("chain B", "if (resi, resn) not in resis:\n\tresis.append((resi, resn))")
start_position = 345
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain B+D+F & resi %s & resn %s"%(r,aa), "resi='%s_'"%(str(i+start_position)))
cmd.iterate("chain B", "print(resi, resn)")

cmd.save("../data/H1_renumbered.pdb")
