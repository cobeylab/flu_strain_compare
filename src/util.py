# Usage:
# cmd.fetch('4we8')
# aa = [i.resn for i in cmd.get_model("chain A" + " and n. ca").atom]
# conv1 = convert_AA_3to1(aa)
# conv3 = convert_AA_1to3(conv1)
# print(conv1, conv3)

from pymol import cmd

# Prepare constants
AA1 = list("ACDEFGHIKLMNPQRSTVWY")
AA3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
AA1to3 = dict(zip(AA1, AA3))
AA3to1 = dict(zip(AA3, AA1))


# For conservative mutations
positive = {'R', 'K', 'H'}
negative = {'D', 'E'}
aromatic = {'F', 'W', 'Y'}
aliphatic = {'A', 'V', 'I', 'L'}
hydroxyl = {'S', 'T'}
amines = {'N', 'Q'}
g = {'G'}
c = {'C'}
p = {'P'}
m = {'M'}


AA_GROUPS = [positive, negative, aromatic, aliphatic, hydroxyl, amines, g, c, p, m]


# Convert three-letter amino acid sequence to one-letter.
def convert_AA_3to1(aa):
    return [AA3to1[i] for i in aa]

# Convert one-letter amino acid sequence to three-letter.
def convert_AA_1to3(aa):
    return [AA1to3[i] for i in aa]

def conservative(prev, curr):
    prev_index = [(prev in s) for s in AA_GROUPS].index(True)
    curr_index = [(curr in s) for s in AA_GROUPS].index(True)
    return prev_index == curr_index

# From https://doi.org/10.3389/fmolb.2015.00056
# A script to highlight hydrophobicity and charge on protein surfaces
# DHS065 Hagemans et al YRB script
# created by Dominique Hagemans and Ianthe A.E.M. van Belzen, July 2015
# Rudiger group CPC, Utrecht University

# yellow: C, CH, CH2, CH3 groups that are not bound to N or O groups.
# red: negatively charged atoms
# blue: positively charged atoms
# grey: backbone, polar groups and remaining atoms

# usage:
# (1) save this script as file with a .py extension at desirable location
# (2) open your structure in pymol and change to surface view.
# (3) run this script: File --> Run --> DHS065 Hagemans et al YRB script
# (4) give the command to colour all structures or a specific structure:
# (4 a) "yrb" to colour all structures
# (4 b) "yrb 'designation'" to colour only that specific structure

def yrb(selection='all'):

	cmd.remove("hydro")
	cmd.set_color('yellow',[0.950,0.78,0.0])
	cmd.set_color('grey',[0.95,0.95,0.95])
	cmd.set_color('red',[1.0,0.4,0.4])
	cmd.set_color('blue',[0.2,0.5,0.8])

	mapping = {}
	mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
	mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
	mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
	mapping['cys'] = [ ('SG', 'grey') ]
	mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
	mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
	mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ]
	mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
	mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
	mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
	mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
	mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
	mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
	mapping['ser'] = [ ('CB,OG', 'grey') ]
	mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
	mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
	mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
	mapping['val'] = [ ('CG1,CG2', 'yellow') ]

	obj_list = cmd.get_names('objects')
	for obj in obj_list:
		if (obj == selection or selection == 'all'):
			cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
			cmd.color('yellow','(n. CB and ' + obj + ')')

			for key in mapping:
				for (atom, color) in mapping[key]:
					cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

