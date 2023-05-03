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
AA_GROUPS = [positive, negative, aromatic, aliphatic, hydroxyl, amines]


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
