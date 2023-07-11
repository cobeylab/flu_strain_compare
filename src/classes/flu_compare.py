import os
from Bio import SeqIO
import re
from os import system
import pandas as pd
from pymol import cmd
import numpy
import json
from json import JSONEncoder


# Internal data directory
DATA_DIR = "data"

class FluMutation:
    def __init__(self,
        pymol_resi,
        label,
        strain1,
        strain2):
        assert (type(pymol_resi) == str), "Position must be a string"
        assert (type(label) == str), "Mutation must be identified as a string"
        assert (type(strain1) == str), "Strain 1 must be a string identifier"
        assert (type(strain2) == str), "Strain 2 must be a string identifier"
        self.pymol_resi = pymol_resi
        self.label = label
        self.strain1 = strain1
        self.strain2 = strain2
    def __str__(self):
        return f"Strain 1: {self.strain1}, Strain 2: {self.strain2}, Mutation: {self.label}"

class FluMutationMultiWay:
    def __init__(self,
        pymol_resi,
        label):
        assert (type(pymol_resi) == str), "Position must be a string"
        assert (type(label) == str), "Mutation must be identified as a string"
        self.pymol_resi = pymol_resi
        self.label = label
    def __str__(self):
        return f"PyMOL residue: {self.pymol_resi}, Mutation: {self.label}"

class FluPngs:
    def __init__(self,
        pymol_resi,
        label,
        percent_conserved):
        assert (type(pymol_resi) == str), "Position must be a string"
        assert (type(label) == str), "Label must be identified as a string"
        self.pymol_resi = pymol_resi
        self.label = label
        self.percent_conserved = percent_conserved

    def __str__(self):
        return f"PNGS: {self.pymol_resi}, Label: {self.label}, Conserved: {self.percent_conserved}"

    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.pymol_resi == other.pymol_resi and self.label == other.label and self.percent_conserved == other.percent_conserved
        return False


class FluSeq:
    def __init__(self,
        lineage,
        query_sequence_file,
        query_sequence_id):
        assert (os.path.exists(query_sequence_file)), f"Query sequence file {query_sequence_file} does not exist"
        assert (type(lineage) == str), "Lineage must be a string."
        assert (type(query_sequence_id) == str), "Query sequence ID must be a string."
        query_seqs = {s.id: s for s in SeqIO.parse(query_sequence_file, "fasta")}
        assert (len(query_seqs) >= 1), "Query sequence file must contain at least 1 sequence."
        self.lineage = lineage
        self.sequence = query_seqs[query_sequence_id]
        self.name = self.sequence.description.split(" | ")[0]
        self.query_sequence_file = query_sequence_file
        self.align_to_reference()

        # Get PNGS sites
        gly = re.compile("N[-]*[A-O,Q-Z][-]*[S,T]")
        self.pngs = [str(m.start() + 1) for m in gly.finditer(str(self.sequence.seq))]

    def align_to_reference(self):
        ref_file = f"{DATA_DIR}/{self.lineage}_ref.fasta"

        TEMP_DIR = f"{DATA_DIR}/tmp"
        if not os.path.exists(TEMP_DIR):
            os.makedirs(TEMP_DIR)

        temp_seqfile = f"{TEMP_DIR}/tmp.fasta"
        temp_alignfile = f"{TEMP_DIR}/aligned.fasta"
        SeqIO.write([self.sequence], temp_seqfile, "fasta")
        command = "mafft --keeplength --add %s %s > %s"%(temp_seqfile, ref_file, temp_alignfile)
        system(command)

        newseq = [s for s in SeqIO.parse(temp_alignfile, "fasta") if s.id == self.sequence.id]
        assert (len(newseq) == 1), f"Alignment {ref_file} not found"
        self.sequence = newseq[0]

        # Clean up temp files
        os.remove(temp_seqfile)
        os.remove(temp_alignfile)

    def __str__(self):
        return f"Name: {self.name}, Lineage: {self.lineage}, Query Sequence File: {self.query_sequence_file}, PNGS: {self.pngs}"

class SequenceComparison:
    def __init__(self, seq1, comparisons, numbering_scheme, reference_mode, filter_sites, reverse_filter_sites):
        assert (type(seq1) == FluSeq), "Seq1 must be a FluSeq object"
        assert (type(comparisons) == list), "comparisons must be a list"
        assert len({s.lineage for s in comparisons} | {seq1.lineage }), "Cannot compare sequences from different lineages"
        self.seq1 = seq1
        self.comparisons = comparisons
        self.lineage = seq1.lineage
        self.numbering_scheme = numbering_scheme
        conversion_file = f"{DATA_DIR}/{self.lineage}_Conversion.csv"
        assert (os.path.exists(conversion_file)), f"Conversion file {conversion_file} does not exist"
        self.conversion_table = pd.read_csv(conversion_file,
            index_col = "ref_one_index")
        self.filter_sites = [str(f) for f in filter_sites]
        self.reverse_filter_sites = [str(f) for f in reverse_filter_sites]

        # Reverse filter takes precedence over filter.
        if len(self.reverse_filter_sites) > 0 and len(self.filter_sites) > 0:
            self.filter_sites = []

        self.mutation_list = self.identify_mutations()
        self.reference_mode = reference_mode
        if self.reference_mode:
            self.gly_del = self.identify_PNGS_changes("deletions")
            self.gly_add = self.identify_PNGS_changes("additions")
            self.gly_share = self.identify_PNGS_changes("shared")
        else:
            self.gly_no_reference = self.identify_PNGS_no_reference()

    def convert_numbering(self,
        position):
        return self.conversion_table.loc[position + 1, self.numbering_scheme]

    def identify_mutations(self):
        # Put sequences in one list
        sequences = [self.seq1.sequence.seq]
        comparisons = [x.sequence.seq for x in self.comparisons]
        sequences.extend(comparisons)

        filter_on = len(self.filter_sites) != 0
        rev_filter_on = len(self.reverse_filter_sites) != 0

        mutations_out = []
        for i, sites in enumerate(zip(*sequences)):
            if not len(set(sites)) == 1:
                p = self.convert_numbering(i)

                if (filter_on and not p in self.filter_sites) or (rev_filter_on and p in self.reverse_filter_sites) or (not filter_on and not_rev_filter_on):
                    mutations_out.append(
                        FluMutationMultiWay(pymol_resi = str(i+1),
                            label = "".join(list(sites) + [str(p)])
                            )
                    )
        return mutations_out

    def identify_PNGS_changes(self, change_type):
        if change_type == "deletions":
            comparison_set = set()
            for comp in self.comparisons:
                comparison_set = comparison_set | set(self.seq1.pngs) - set(comp.pngs)
        elif change_type == "additions":
            comparison_set = set()
            for comp in self.comparisons:
                comparison_set = comparison_set | set(comp.pngs) - set(self.seq1.pngs).intersection(set(comp.pngs))
        elif change_type == "shared":
            comparison_set = set(self.seq1.pngs)
            for comp in self.comparisons:
                comparison_set = comparison_set.intersection(comp.pngs)

        # filter out sites named in filter_sites setting
        comparison_set = comparison_set - set(self.filter_sites)

        if len(self.reverse_filter_sites):
            comparison_set = set(self.reverse_filter_sites) - comparison_set


        pngs_out = []
        for pymol_position in comparison_set:
            conversion_index = int(pymol_position.replace("_","")) - 1
            label = "PNGS%s"%self.convert_numbering(conversion_index)
            pngs_out.append(FluPngs(pymol_resi = pymol_position,
                label = label, percent_conserved = None))
        return pngs_out

    def identify_PNGS_no_reference(self):
        all_comparisons = self.comparisons + [self.seq1]
        all_comparisons_pngs = [set(c.pngs) for c in all_comparisons]
        return compare_seq_no_reference(all_comparisons_pngs, self.convert_numbering, set(self.filter_sites), set(self.reverse_filter_sites))

# Encoder for SequenceComparison serialization
class SequenceComparisonEncoder(JSONEncoder):
        def default(self, o):
            if isinstance(o, numpy.ndarray):
                return o.tolist()
            if isinstance(o, SequenceComparison):
                # Don't serialize conversion table
                sc = o.__dict__
                sc.pop("conversion_table")
                return sc
            return o.__dict__

def compare_seq_no_reference(comparisons, convert, filtered=set(), rev_filtered=set()):
    comparison_set = set()
    for comp in comparisons:
        comparison_set = comparison_set | comp

    # filter out sites named in filter_sites setting
    comparison_set = comparison_set - filtered

    if len(rev_filtered) > 0:
        comparison_set = rev_filtered - comparison_set


    pngs_out = []
    num_comparisons = len(comparison_set)
    for pymol_position in comparison_set:
        conversion_index = int(pymol_position.replace("_","")) - 1
        label = "PNGS%s"%convert(conversion_index)
        percent_conserved = sum([pymol_position in comp for comp in comparisons])/num_comparisons
        pngs_out.append(FluPngs(pymol_resi = pymol_position,
            label = label, percent_conserved = percent_conserved))

    return pngs_out


# Render a SequenceComparison
def make_figure(sc):
    names = [sc.seq1.name] + [s.name for s in sc.comparisons]
    mutations = [m.pymol_resi for m in sc.mutation_list if m.pymol_resi != "-"]
    cmd.reinitialize()
    cmd.set('ray_trace_mode', 0)

    # Label parameters
    cmd.set('label_position', (0, 0, 20))
    cmd.set('label_size', -3)
    cmd.set('label_color', 'black')
    cmd.load(f"{DATA_DIR}/{sc.lineage}_pngs.pse")

    # Color mutations
    if len(mutations) > 0:
        cmd.select('mutations', '(resi %s)'%'+'.join([i for i in mutations]))
        cmd.color('yellow', 'mutations')

    # Color glycosylations
    if sc.reference_mode:
        color_pngs(sc.gly_del, "glycan_deletions", "red")
        color_pngs(sc.gly_add, "glycan_additions", "green")
        color_pngs(sc.gly_share, "glycans_shared", "blue")
    else:
        color_pngs_no_reference(sc.gly_no_reference, "glycans_no_reference", "red")

    # Name file with the first 3 strains.
    base_filename = "-".join([n.replace("/", "_") for n in names[:3]])
    names_len = len(names)
    if names_len > 3:
        base_filename += f"-{names_len - 4}_others"

    # Add labels
    label_resi(sc.mutation_list)
    if sc.reference_mode:
        label_resi(sc.gly_del)
        label_resi(sc.gly_add)
        label_resi(sc.gly_share)

        x_pos = -60
        z_pos = 0
        y_start = 0
        y_offset = -5
        create_label(x_pos, y_start, z_pos, "Mutations", "mutlabel", "yellow")
        create_label(x_pos, y_start + y_offset, z_pos, "PNGS added", "glyaddlabel", "green")
        create_label(x_pos, y_start + 2*y_offset , z_pos, "PNGS deleted", "glydellabel", "red")
        create_label(x_pos, y_start + 3*y_offset, z_pos, "PNGS shared", "glysharelabel", "blue")
    else:
        label_resi(sc.gly_no_reference)

    cmd.hide("everything", "extra_glycans")

    create_label(0, 110, 0, " vs. ".join(names), "strains", "white", label_size=-5)
    return base_filename

def color_pngs(glylist, name, color):
    if len(glylist) > 0:
        PNGS_names = ["PNGS%s"%g.pymol_resi for g in glylist]
        PNGS_names_final = set(PNGS_names).intersection(set(cmd.get_names(type="selections")))
        # Need to add warning here if it finds a PNGS site that isn't in the structure
        if len(PNGS_names_final) > 0:
            cmd.select(name, ' | '.join(PNGS_names_final))
            cmd.show("sticks", name)
            cmd.color(color, name)

def color_pngs_no_reference(glylist, name, color):
    if len(glylist) > 0:
        PNGS_names = ["PNGS%s"%g.pymol_resi for g in glylist]
        PNGS_names_final = set(PNGS_names).intersection(set(cmd.get_names(type="selections")))
        # Need to add warning here if it finds a PNGS site that isn't in the structure
        if len(PNGS_names_final) > 0:
            cmd.select(name, ' | '.join(PNGS_names_final))
            cmd.show("sticks", name)

        for g in glylist:
            try:
                cmd.select(name + g.pymol_resi, "PNGS%s"%g.pymol_resi)
                cmd.set_color("color_" + g.pymol_resi, [1 - g.percent_conserved, 0.0, 0.0])
                cmd.color("color_" + g.pymol_resi, name + g.pymol_resi)
            except Exception as e:
                print(e)

def create_label(x, y, z, label_text, label_name, label_color, label_size=-4):
    cmd.pseudoatom(label_name)
    cmd.label(selection = label_name, expression = f"'{label_text}'")
    cmd.hide("wire", selection = label_name)
    cmd.set("label_color", selection = label_name, value = label_color)
    cmd.set("label_position", selection = label_name, value = [x,y,z])
    cmd.set("label_size", selection = label_name, value = label_size)

def label_resi(resilist):
    for m in resilist:
        label = m.label
        resi = m.pymol_resi
        if resi != "-":
            cmd.select(label, 'n. CA and i. ' + resi)
            cmd.label(selection = label, expression = f"'{label}'")

def make_comparison_object(parameters):
    seq_file = parameters["seq_file"]
    reference_id = parameters["reference_id"]
    comparison_ids = parameters["comparison_ids"]
    seq_lineage = parameters["seq_lineage"]
    numbering_scheme = parameters["numbering_scheme"]
    reference_mode = parameters["reference_mode"]
    filter_sites = parameters.get("filter_sites", [])
    reverse_filter_sites = parameters.get("reverse_filter_sites", [])

    s1 = FluSeq(
        lineage = seq_lineage,
        query_sequence_file = seq_file,
        query_sequence_id = reference_id,
        )
    comparisons = []
    for c_id in comparison_ids:
        seq = FluSeq(
            lineage = seq_lineage,
            query_sequence_file = seq_file,
            query_sequence_id = c_id,
            )

        comparisons.append(seq)

    comparison = SequenceComparison(seq1 = s1,
        comparisons = comparisons,
        numbering_scheme = numbering_scheme,
        reference_mode = reference_mode,
        filter_sites = filter_sites,
        reverse_filter_sites = reverse_filter_sites
        )

    return comparison

