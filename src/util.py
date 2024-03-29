# Usage:
# cmd.fetch('4we8')
# aa = [i.resn for i in cmd.get_model("chain A" + " and n. ca").atom]
# conv1 = convert_AA_3to1(aa)
# conv3 = convert_AA_1to3(conv1)
# print(conv1, conv3)

from pymol import cmd
from scipy.stats import entropy
from math import log

# Prepare constants
AA1 = list("ACDEFGHIKLMNPQRSTVWY")
AA3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
AA1to3 = dict(zip(AA1, AA3))
AA3to1 = dict(zip(AA3, AA1))


# For conservative mutations
positive = ['R', 'K', 'H']
negative = ['D', 'E']
aromatic = ['F', 'W', 'Y']
aliphatic = ['A', 'V', 'I', 'L']
hydroxyl = ['S', 'T']
amines = ['N', 'Q']
g = ['G']
c = ['C']
p = ['P']
m = ['M']


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

# Convert a set of mutations to conservative. (Replace all amino acids with the first element of their AA_GROUP.)
def to_conservative(sites):
    new_sites = []
    for s in sites:
        index = [(s in group) for group in AA_GROUPS].index(True)
        group = AA_GROUPS[index]
        new_sites.append(group[0])
    return new_sites

# Shannon index normalized to the range 0.0-1.0.
def shannon(events):
    num_events = len(events)
    unique_events = set(events)
    probs = [events.count(e)/num_events for e in unique_events]
    return entropy(probs)/log(num_events)

# Richness normalized to the range 0.0-1.0.
def richness(events):
    return len(set(events))/len(events)

# Gini-Simpson Index
def gini_simpson(events):
    gsi = 1 - sum([pow(proportional_abundance(e, events), 2) for e in set(events)])
    return gsi

# Proportional abundance (probability of selecting an element).
def proportional_abundance(element, events):
    return events.count(element)/len(events)
