from pymol import cmd
cmd.reinitialize()
cmd.fetch('4we8') # A/Victoria/361/2011 sequence for reference

# These next lines of code expand the monomer to a trimer
cmd.symexp('sym','4we8','4we8','3')
cmd.delete('sym11000000')
cmd.delete('sym04000000')
cmd.center('4we8, sym01000000, sym02000000')

# Basic visual settings
cmd.show('surface')
cmd.remove('solvent')
cmd.color('gray70', 'all')
cmd.rotate("y", "80")
cmd.rotate("z", "-90")
cmd.rotate("x", "40")
cmd.rotate("y", "25")
cmd.rotate("x", "10")
cmd.rotate("z", "20")

# Rremove glycans
# In the structure 4we8, the glycans can be identified as NAG, MAN, and BMA
cmd.remove("resn NAG+MAN+BMA")

# The loop below renumbers all the residues so that they match
# Standard H3N2 numbering with position 1 beginning after the leader
# peptide.

resis = []
cmd.iterate("4we8", "if (resi, resn) not in resis:\n\tresis.append((resi, resn))")
offset = 16
for i, (r, aa) in enumerate(resis):
	if aa != "NAG":
		new_resi = str(int(r)+offset)
		cmd.alter("resi %s and resn %s"%(r,aa), "resi='%s_'"%(new_resi))
cmd.iterate("4we8", "print(resi, resn)")
cmd.save("data/H3_renumbered.pse")
