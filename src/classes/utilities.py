from pymol import cmd
import sys
sys.path.append("/app/src/classes")
from flu_compare import flu_seq,seq_compare

def color_pngs(glylist, name, color):
    if len(glylist) > 0:
        cmd.select(name, '(resi %s)'%'+'.join([g.pymol_resi for g in glylist]))
        cmd.color(color, name)
def create_label(x, y, z, label_text, label_name, label_color):
    cmd.pseudoatom(label_name, pos=[x, y, z])
    global label_temp
    label_temp = label_text
    cmd.label(selection = label_name, expression = "label_temp")
    cmd.hide("wire", selection = label_name)
    cmd.set("label_color", selection = label_name, value = label_color)
def label_resi(resilist):
    for m in resilist:
        label = m.label
        resi = m.pymol_resi
        if resi != "-":
            cmd.select(label, 'n. CA and i. ' + resi)
            global label_name
            label_name = str(label)
            cmd.label(selection = label, expression = "label_name")
def make_comparison_object(parameters):
    seq_file = "/app/data/" + parameters["seq_file"]
    q1_id = parameters["q1_id"]
    q1_name = parameters["q1_name"]
    q2_id = parameters["q2_id"]
    q2_name = parameters["q2_name"]
    seq_lineage = parameters["seq_lineage"]
    numbering_scheme = parameters["numbering_scheme"]
    s1 = flu_seq(name = q1_name,
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q1_id
        )
    s2 = flu_seq(name = q2_name,
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = q2_id
        )
    comparison = seq_compare(seq1 = s1, seq2 = s2, numbering_scheme=numbering_scheme)
    mutation_list = comparison.identify_mutations()
    mutations = [m.pymol_resi for m in mutation_list if m.pymol_resi != "-"]
    gly_del = comparison.identify_PNGS_changes("deletions")
    gly_add = comparison.identify_PNGS_changes("additions")
    gly_share = comparison.identify_PNGS_changes("shared")
    return (mutations, gly_del, gly_add, gly_share, mutation_list, q1_name, q2_name, seq_lineage)
