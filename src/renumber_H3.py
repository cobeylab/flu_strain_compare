from pymol import cmd

remove_glycans = False

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
if (remove_glycans):	
	cmd.remove("resn NAG+MAN+BMA")
	suffix = "no_pngs"
else:
	suffix = "with_pngs"

# The loop below renumbers all the residues so that they match
# Standard H3N2 numbering with position 1 beginning after the leader
# peptide.

# In newer versions of pymol, you need to explicity specify which variables 
# will get passed into cmd
myspace = {'resis': []}
resis = myspace['resis']

cmd.iterate("4we8", "if (resi, resn) not in resis:\n\tresis.append((resi, resn))", space = myspace)

offset = 16
for i, (r, aa) in enumerate(resis):
	if aa not in ["NAG", "MAN", "BMA"]:
		new_resi = str(int(r)+offset)
		cmd.alter("resi %s and resn %s"%(r,aa), "resi='%s_'"%(new_resi))
cmd.iterate("4we8", "print(resi, resn)")
cmd.save("H3_renumbered_%s.pse"%(suffix))