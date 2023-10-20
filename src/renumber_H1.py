from pymol import cmd

remove_glycans = False

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

if (remove_glycans):	
	cmd.remove("resi 400+401+402+403+407")
	suffix = "no_pngs"
else:
	suffix = "with_pngs"

# Removing glycans
# In this case, it was easier to remove residues that represented glycans and not actual amino acids
# I figured this out manually when I first played around with this structure
cmd.remove("chain G+H")

# The following loops renumber the residues so that position 1 indicates the first
# residue after the leader sequence.

myspace = {'resis': [], 'excluded_resi': ['400', '401', '402', '403', '407']}
resis = myspace['resis']
cmd.iterate("chain A", "if (resi, resn) not in resis and resi not in excluded_resi:\n\tresis.append((resi, resn))", space = myspace)
start_position = 14
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain A+C+E & resi %s & resn %s"%(r,aa), "resi='%s_'"%(new_resi))

cmd.iterate("chain A", "print(resi, resn)")


myspace = {'resis': [], 'excluded_resi': ['400', '401', '402', '403', '407']}
resis = myspace['resis']
cmd.iterate("chain B", "if (resi, resn) not in resis and resi not in excluded_resi:\n\tresis.append((resi, resn))", space = myspace)
start_position = 345
for i, (r, aa) in enumerate(resis):
	new_resi = str(i+start_position)
	cmd.alter("chain B+D+F & resi %s & resn %s"%(r,aa), "resi='%s_'"%(str(i+start_position)))
cmd.iterate("chain B", "print(resi, resn)")

cmd.save("H1_renumbered_%s.pse"%(suffix))
