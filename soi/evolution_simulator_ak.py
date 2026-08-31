#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Rearrangement Evolution Simulator (AK.py output format)
==================================================================

Simulates realistic chromosome rearrangements along a phylogenetic tree,
producing output compatible with AK.py (Ancestral Karyotype Reconstruction).

Events: Inversion, RT, NCF, EEJ, Fission, WGD, Fractionation,
        Tandem duplication, Dispersed duplication, Gene gain, Gene loss,
        Unidirectional translocation,
        Segmental deletion, Segmental duplication, Chromothripsis.

Usage:
    soi sim -t tree.nwk -o output_dir [options]
    soi sim --num-species 6 -o output_dir [options]   # auto-generate tree
"""

from __future__ import print_function
import sys
import os
import copy
import math
import random
import string
import argparse
import re
from collections import defaultdict, Counter

from ete3 import Tree
import networkx as nx

from .ak_dotplot import draw_ancestor_dotplots, draw_branch_dotplots, draw_sister_dotplots
from .tree import convert_newick


# ============================================================================
# Helper functions
# ============================================================================

Gene = tuple  # (gene_id: str, orientation: str '+'/'-')
Karyotype = dict  # {chrom_id: [Gene, ...]}


def flip_orient(o):
    return "+" if o == "-" else "-"


def reverse_segment(genes):
    return [(g, flip_orient(o)) for g, o in reversed(genes)]


def get_chr_labels(n):
    return [str(i + 1) for i in range(n)]


def poisson_sample(rng, lam):
    if lam <= 0:
        return 0
    if lam < 30:
        L = math.exp(-lam)
        k = 0
        p = 1.0
        while True:
            k += 1
            p *= rng.random()
            if p <= L:
                return k - 1
    else:
        val = rng.gauss(lam, math.sqrt(lam))
        return max(0, int(round(val)))


def sample_wgd_factor(rng):
    """P(2) = 1/2, P(3) = 1/4, P(4) = 1/8, ..."""
    r = rng.random()
    n = 2
    cum = 0.5
    while r > cum:
        n += 1
        cum += 0.5 ** (n - 1)
    return n


def generate_tree(num_species, rng):
    """Generate a random balanced binary tree with branch lengths."""
    names = ["Sp_{}".format(i) for i in range(1, num_species + 1)]
    rng.shuffle(names)
    nodes = names[:]
    while len(nodes) > 1:
        rng.shuffle(nodes)
        next_nodes = []
        i = 0
        while i < len(nodes):
            if i + 1 < len(nodes):
                bl1 = round(rng.uniform(0.05, 0.30), 3)
                bl2 = round(rng.uniform(0.05, 0.30), 3)
                pair = "({}:{},{}:{})".format(nodes[i], bl1, nodes[i+1], bl2)
                next_nodes.append(pair)
            else:
                next_nodes.append(nodes[i])
            i += 2
        nodes = next_nodes
    nw = nodes[0] + ";"
    # Parse and add root branch length
    tree = Tree(nw, format=1)
    return tree


def parse_tree(tree_file, wgd_rate=0.0):
    """Parse tree file, extracting [p=N] ploidy annotations from any node.
    Tree annotations always take precedence; wgd_rate only controls
    additional rate-based WGD on unannotated nodes."""
    with open(tree_file) as f:
        orig_nw = f.read().strip()

    # Convert [p=N] annotations to NHX format so ete3 can parse them
    # as node features (works for leaves, named internals, unnamed internals)
    converted_nw = convert_newick(orig_nw)

    tree = None
    for fmt in [1, 0, 2, 3, 8]:
        try:
            tree = Tree(converted_nw, format=fmt)
            break
        except Exception:
            continue
    if tree is None:
        raise ValueError("Cannot parse tree file: " + tree_file)

    existing_names = set(node.name for node in tree.traverse() if node.name)
    counter = 0
    for node in tree.traverse("preorder"):
        if not node.name or node.name.strip() == "":
            while True:
                candidate = "N{}".format(counter)
                counter += 1
                if candidate not in existing_names:
                    node.name = candidate
                    existing_names.add(candidate)
                    break

    # Build ploidy_map from ete3 node features (NHX:p=N)
    ploidy_map = {}
    for node in tree.traverse():
        if hasattr(node, 'p'):
            try:
                ploidy = int(node.p)
                if ploidy > 1:
                    ploidy_map[node.name] = ploidy
            except (ValueError, TypeError):
                pass

    return tree, ploidy_map, orig_nw


def number_tree_nodes(tree):
    """Number internal tree nodes as N0, N1, ..."""
    i = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        if node.name and re.match(r'^N\d+$', node.name):
            continue
        node.name = "N{}".format(i)
        i += 1
    internal_nodes = [n for n in tree.traverse("postorder") if not n.is_leaf()]
    for idx, node in enumerate(internal_nodes):
        node.name = "N{}".format(idx)


# ============================================================================
# GeneTracker -- dot-notation naming for homologs
# ============================================================================

class GeneTracker:
    def __init__(self):
        self._anc_counter = 0
        self._gain_counter = 0
        self._homolog_seq = {}     # anc_gene -> next suffix index
        self.ancestor_genes = []
        self.gene_ancestor = {}
        self.gene_species = {}
        self.gene_wgd_copies = {}
        self.gene_is_gain = {}
        self.gene_dup_of = {}
        self.gene_dup_type = {}
        self.gene_parent = {}       # maps gene_id -> immediate parent gene_id
        self.gain_genes = []

    def _next_homolog_id(self, anc_gene):
        idx = self._homolog_seq.get(anc_gene, 1)
        self._homolog_seq[anc_gene] = idx + 1
        return "{}.{}".format(anc_gene, idx)

    def new_ancestral_gene(self):
        self._anc_counter += 1
        gid = "g{}".format(self._anc_counter)
        self.gene_ancestor[gid] = gid
        self.ancestor_genes.append(gid)
        self.gene_wgd_copies[gid] = {}
        self.gene_is_gain[gid] = False
        return gid

    def new_wgd_copy(self, original_gid, wgd_node, copy_index, species=None):
        anc = self.gene_ancestor.get(original_gid, original_gid)
        gid = self._next_homolog_id(anc)
        self.gene_ancestor[gid] = anc
        self.gene_parent[gid] = original_gid
        if species is not None:
            self.gene_species[gid] = species
        wgd = dict(self.gene_wgd_copies.get(original_gid, {}))
        wgd[wgd_node] = copy_index
        self.gene_wgd_copies[gid] = wgd
        self.gene_is_gain[gid] = False
        return gid

    def new_duplicate_gene(self, original_gid, dup_type, species=None):
        anc = self.gene_ancestor.get(original_gid, original_gid)
        gid = self._next_homolog_id(anc)
        self.gene_ancestor[gid] = anc
        self.gene_parent[gid] = original_gid
        if species is not None:
            self.gene_species[gid] = species
        self.gene_wgd_copies[gid] = dict(self.gene_wgd_copies.get(original_gid, {}))
        self.gene_is_gain[gid] = False
        self.gene_dup_of[gid] = original_gid
        self.gene_dup_type[gid] = dup_type
        return gid

    def new_gain_gene(self, species):
        self._gain_counter += 1
        gid = "gn{}".format(self._gain_counter)
        self.gene_ancestor[gid] = gid
        self.gene_species[gid] = species
        self.gene_wgd_copies[gid] = {}
        self.gene_is_gain[gid] = True
        self.gain_genes.append(gid)
        return gid

    def get_ancestor(self, gid):
        return self.gene_ancestor.get(gid, gid)

    def set_species(self, gid, species):
        self.gene_species[gid] = species

    def total_genes(self):
        total = self._anc_counter + self._gain_counter
        total += sum(v - 1 for v in self._homolog_seq.values())
        return total


# ============================================================================
# EvolutionSimulator
# ============================================================================

class EvolutionSimulator:
    def __init__(self, seed=42, num_chroms=8, min_genes=200, max_genes=1000,
                 inv_rate=10.0, rt_rate=3.0, ncf_rate=1.5, eej_rate=2.0,
                 fission_rate=0.02, wgd_rate=0.4,
                 gene_gain_rate=0.5, gene_loss_rate=0.5, tandem_dup_rate=2.0,
                 dispersed_dup_rate=0.3, unidir_trans_rate=1.0,
                 frac_rate=5.0, scale=1.0,
                 seg_del_rate=0.5, seg_dup_rate=0.3, chromothripsis_rate=0.05,
                 inv_mode="random", rt_mode="random", eej_mode="random"):
        self.seed = seed
        self.rng = random.Random(seed)
        self.num_chroms = num_chroms
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.inv_rate = inv_rate * scale
        self.rt_rate = rt_rate * scale
        self.ncf_rate = ncf_rate * scale
        self.eej_rate = eej_rate * scale
        self.fission_rate = fission_rate * scale
        self.wgd_rate = wgd_rate * scale
        self.gene_gain_rate = gene_gain_rate * scale
        self.gene_loss_rate = gene_loss_rate * scale
        self.tandem_dup_rate = tandem_dup_rate * scale
        self.dispersed_dup_rate = dispersed_dup_rate * scale
        self.unidir_trans_rate = unidir_trans_rate * scale
        self.frac_rate = frac_rate * scale
        self.seg_del_rate = seg_del_rate * scale
        self.seg_dup_rate = seg_dup_rate * scale
        self.chromothripsis_rate = chromothripsis_rate * scale
        self.rt_mode = rt_mode
        self.eej_mode = eej_mode
        self.inv_mode = inv_mode

        self.tracker = GeneTracker()
        self.leaf_karyotypes = {}
        self.all_node_karyotypes = {}
        self.events = []
        self.wgd_map = {}
        self.centromeres = {}
        self.wgd_map = {}
        self.centromeres = {}

    def run(self, tree, ploidy_map):
        print("Initializing ancestor with {} chromosomes ({}-{} genes each)...".format(
            self.num_chroms, self.min_genes, self.max_genes))
        root_karyo, root_centros = self._init_ancestor()
        n_leaves = len(tree.get_leaves())
        print("Simulating evolution along tree with {} species...".format(n_leaves))
        self._evolve(tree, root_karyo, root_centros, ploidy_map)
        self._print_summary()
        return self.leaf_karyotypes

    def _init_ancestor(self):
        chr_labels = get_chr_labels(self.num_chroms)
        karyo = {}
        centros = {}
        for cid in chr_labels:
            n = self.rng.randint(self.min_genes, self.max_genes)
            genes = [(self.tracker.new_ancestral_gene(), "+") for _ in range(n)]
            karyo[cid] = genes
            centros[cid] = self.rng.randint(2, n - 2) if n > 4 else n // 2
        return karyo, centros

    def _evolve(self, tree, root_karyo, root_centros, ploidy_map):
        for node in tree.traverse("preorder"):
            if node.is_root():
                self.all_node_karyotypes[node.name] = root_karyo
                self.centromeres[node.name] = root_centros
                continue

            karyo = copy.deepcopy(self.all_node_karyotypes[node.up.name])
            centros = copy.deepcopy(self.centromeres[node.up.name])
            name = node.name
            bl = node.dist if node.dist else 0.0

            # WGD — save pre-WGD snapshot before applying
            if name in ploidy_map and ploidy_map[name] > 1:
                pre_wgd_name = "{}_preWGD".format(name)
                self.all_node_karyotypes[pre_wgd_name] = copy.deepcopy(karyo)
                self.centromeres[pre_wgd_name] = copy.deepcopy(centros)
                factor = ploidy_map[name]
                print("  Applying {}x genome duplication at {} (from tree annotation)".format(
                    factor, name))
                karyo, centros = self._apply_wgd(karyo, centros, name, factor)
            elif self.wgd_rate > 0 and bl > 0:
                n_wgd = poisson_sample(self.rng, self.wgd_rate * bl)
                for _ in range(n_wgd):
                    pre_wgd_name = "{}_preWGD".format(name)
                    if pre_wgd_name not in self.all_node_karyotypes:
                        self.all_node_karyotypes[pre_wgd_name] = copy.deepcopy(karyo)
                        self.centromeres[pre_wgd_name] = copy.deepcopy(centros)
                    factor = sample_wgd_factor(self.rng)
                    print("  Applying {}x genome duplication at {} (rate-based)".format(
                        factor, name))
                    karyo, centros = self._apply_wgd(karyo, centros, name, factor)

            if bl > 0:
                self._apply_rearrangements(karyo, centros, name, bl)

            self.all_node_karyotypes[name] = karyo
            self.centromeres[name] = centros
            if node.is_leaf():
                self.leaf_karyotypes[name] = karyo
                for cid, genes in karyo.items():
                    for gid, orient in genes:
                        self.tracker.set_species(gid, name)

    # ------------------------------------------------------------------
    # WGD
    # ------------------------------------------------------------------
    def _apply_wgd(self, karyo, centros, node_name, factor):
        new_k = {}
        new_centros = {}
        for c, g in karyo.items():
            for copy_i in range(1, factor + 1):
                new_c = "{}_{}".format(c, copy_i)
                new_g = []
                for gid, orient in g:
                    new_gid = self.tracker.new_wgd_copy(gid, node_name, copy_i)
                    new_g.append((new_gid, orient))
                new_k[new_c] = new_g
                new_centros[new_c] = centros.get(c, len(new_g) // 2)
        self.events.append({"node": node_name, "type": "WGD", "factor": factor,
                            "desc": "{}->{} chromosomes".format(len(karyo), len(new_k))})
        # Track max ploidy per species
        cur = self.wgd_map.get(node_name, 1)
        self.wgd_map[node_name] = max(cur, factor)
        return new_k, new_centros

    # ------------------------------------------------------------------
    # Fractionation
    # ------------------------------------------------------------------
    def _apply_fractionation(self, karyo, centros, node_name):
        # Find all multi-copy genes (same ancestor has >=2 copies)
        anc_copies = defaultdict(list)
        for cid in list(karyo.keys()):
            for i, (gid, orient) in enumerate(karyo[cid]):
                anc = self.tracker.get_ancestor(gid)
                anc_copies[anc].append((gid, orient, cid, i))

        # Only keep ancestors with >1 copy
        multi = {anc: copies for anc, copies in anc_copies.items()
                 if len(copies) > 1}
        if not multi:
            return False

        # Randomly pick one ancestor and remove one random copy
        anc = self.rng.choice(list(multi.keys()))
        copies = multi[anc]
        victim = self.rng.choice(copies)
        lost_gid, lost_orient, lost_cid, lost_idx = victim

        # Remove from karyotype
        genes = karyo[lost_cid]
        retained = [g for i, g in enumerate(genes) if i != lost_idx]
        # Update centromere
        if lost_cid in centros:
            c = centros[lost_cid]
            if c > lost_idx:
                centros[lost_cid] = c - 1
        if retained:
            karyo[lost_cid] = retained
        else:
            del karyo[lost_cid]
            centros.pop(lost_cid, None)

        # Remove from tracker
        self.tracker.gene_ancestor.pop(lost_gid, None)
        self.tracker.gene_species.pop(lost_gid, None)
        self.tracker.gene_wgd_copies.pop(lost_gid, None)

        self.events.append({"node": node_name, "type": "fractionation",
                            "gene": lost_gid, "chrom": lost_cid,
                            "anc": anc, "len": 1})
        return True

    # ------------------------------------------------------------------
    # Rearrangement dispatch
    # ------------------------------------------------------------------
    def _apply_rearrangements(self, karyo, centros, node_name, branch_length):
        n_chrom = len(karyo)
        chrom_factor = n_chrom / max(1, self.num_chroms)
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length * chrom_factor)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length * chrom_factor)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_loss = poisson_sample(self.rng, self.gene_loss_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)
        n_frac = poisson_sample(self.rng, self.frac_rate * branch_length)
        n_sdel = poisson_sample(self.rng, self.seg_del_rate * branch_length)
        n_sdup = poisson_sample(self.rng, self.seg_dup_rate * branch_length)
        n_cht = poisson_sample(self.rng, self.chromothripsis_rate * branch_length)

        event_list = (
            ["inv"] * n_inv + ["rt"] * n_rt + ["ncf"] * n_ncf +
            ["eej"] * n_eej + ["fis"] * n_fis + ["gain"] * n_gain +
            ["loss"] * n_loss + ["tdup"] * n_tdup + ["ddup"] * n_ddup +
            ["utrans"] * n_utrans + ["frac"] * n_frac +
            ["segdel"] * n_sdel + ["segdup"] * n_sdup + ["cht"] * n_cht
        )
        if not event_list:
            return
        self.rng.shuffle(event_list)

        applied = defaultdict(int)
        for e in event_list:
            if e == "inv" and self._apply_inversion(karyo, node_name, centros):
                applied["inv"] += 1
            elif e == "rt" and self._apply_rt(karyo, node_name, centros):
                applied["rt"] += 1
            elif e == "ncf" and self._apply_ncf(karyo, node_name, centros):
                applied["ncf"] += 1
            elif e == "eej" and self._apply_eej(karyo, node_name, centros):
                applied["eej"] += 1
            elif e == "fis" and self._apply_fission(karyo, node_name, centros):
                applied["fis"] += 1
            elif e == "gain" and self._apply_gene_gain(karyo, node_name, centros):
                applied["gain"] += 1
            elif e == "loss" and self._apply_gene_loss(karyo, node_name, centros):
                applied["loss"] += 1
            elif e == "tdup" and self._apply_tandem_dup(karyo, node_name, centros):
                applied["tdup"] += 1
            elif e == "ddup" and self._apply_dispersed_dup(karyo, node_name, centros):
                applied["ddup"] += 1
            elif e == "utrans" and self._apply_unidir_trans(karyo, node_name, centros):
                applied["utrans"] += 1
            elif e == "frac" and self._apply_fractionation(karyo, centros, node_name):
                applied["frac"] += 1
            elif e == "segdel" and self._apply_seg_deletion(karyo, node_name, centros):
                applied["segdel"] += 1
            elif e == "segdup" and self._apply_seg_duplication(karyo, node_name, centros):
                applied["segdup"] += 1
            elif e == "cht" and self._apply_chromothripsis(karyo, node_name, centros):
                applied["cht"] += 1

        self.events.append({"node": node_name, "type": "rearrangements",
                            "branch_length": branch_length,
                            "sampled": {"inv": n_inv, "rt": n_rt, "ncf": n_ncf,
                                        "eej": n_eej, "fis": n_fis, "gain": n_gain,
                                        "loss": n_loss, "tdup": n_tdup, "ddup": n_ddup,
                                        "utrans": n_utrans, "frac": n_frac,
                                        "segdel": n_sdel, "segdup": n_sdup, "cht": n_cht},
                            "applied": dict(applied)})

    def _apply_inversion(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 2]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        n = len(genes)
        mode = self.inv_mode
        if mode == "telomere":
            # 端粒倒位：p1=0 或 p2=n
            if self.rng.random() < 0.5:
                p1 = 0
                p2 = self.rng.randint(2, n)
            else:
                p1 = self.rng.randint(0, n - 2)
                p2 = n
        elif mode == "internal":
            # 内部倒位：p1>0 且 p2<n
            p1 = self.rng.randint(1, n - 2)
            p2 = self.rng.randint(p1 + 1, n - 1)
        else:
            p1, p2 = sorted(self.rng.sample(range(n), 2))
        if p1 == p2:
            return False
        karyo[cid] = genes[:p1] + reverse_segment(genes[p1:p2]) + genes[p2:]
        self.events.append({"node": node_name, "type": "inversion",
                            "chrom": cid, "pos": "{}-{}".format(p1, p2),
                            "mode": mode})
        return True

    def _apply_rt(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 4]
        if len(cids) < 2:
            return False
        c1, c2 = self.rng.sample(cids, 2)
        g1, g2 = karyo[c1], karyo[c2]
        b1 = self.rng.randint(2, len(g1) - 2)
        b2 = self.rng.randint(2, len(g2) - 2)
        aL, aR = g1[:b1], g1[b1:]
        bL, bR = g2[:b2], g2[b2:]
        mode = self.rt_mode
        if mode == "random":
            mode = self.rng.choice(["tailswap", "alternate"])
        if mode == "alternate":
            karyo[c1] = aL + reverse_segment(bL)
            karyo[c2] = reverse_segment(aR) + bR
        else:
            karyo[c1] = aL + bR
            karyo[c2] = bL + aR
        self.events.append({"node": node_name, "type": "RT",
                            "chroms": [c1, c2],
                            "pos": "{}:{}|{}:{}".format(c1, b1, c2, b2),
                            "mode": mode})
        return True

    def _apply_ncf(self, karyo, node_name, centros=None):
        if len(karyo) < 2:
            return False
        keys = list(karyo.keys())
        recipient_cid, donor_cid = self.rng.sample(keys, 2)
        recipient = karyo[recipient_cid]
        donor = list(karyo[donor_cid])
        if len(recipient) < 2:
            return False
        insert_pos = self.rng.randint(1, len(recipient) - 1)
        if self.rng.choice([True, False]):
            donor = reverse_segment(donor)
        karyo[recipient_cid] = recipient[:insert_pos] + donor + recipient[insert_pos:]
        del karyo[donor_cid]
        if centros:
            if recipient_cid in centros:
                c = centros[recipient_cid]
                if c > insert_pos:
                    centros[recipient_cid] = c + len(donor)
            centros.pop(donor_cid, None)
        self.events.append({"node": node_name, "type": "NCF",
                            "chroms": [recipient_cid, donor_cid],
                            "pos": "insert {} into {}:{}".format(donor_cid, recipient_cid, insert_pos)})
        return True

    def _apply_eej(self, karyo, node_name, centros=None):
        if len(karyo) < 2:
            return False
        c1, c2 = self.rng.sample(list(karyo.keys()), 2)
        g1, g2 = list(karyo[c1]), list(karyo[c2])
        e1 = self.rng.choice(["H", "T"])
        e2 = self.rng.choice(["H", "T"])
        if e1 == "H":
            g1 = reverse_segment(g1)
        if e2 == "T":
            g2 = reverse_segment(g2)
        karyo[c1] = g1 + g2
        del karyo[c2]
        if centros:
            if c1 in centros:
                c = centros[c1]
                if e1 == "H":
                    c = len(g1) - c
                centros[c1] = c
            centros.pop(c2, None)
        self.events.append({"node": node_name, "type": "EEJ",
                            "chroms": [c1, c2], "mode": "{}{}".format(e1, e2)})
        return True

    def _apply_fission(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 4]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        # Centromere constraint: fission must occur at centromere
        if centros and cid in centros:
            pos = centros[cid]
            if pos < 2 or pos > len(genes) - 2:
                return False
        else:
            pos = self.rng.randint(2, len(genes) - 2)
        new_id = cid + "f"
        while new_id in karyo:
            new_id += "f"
        karyo[new_id] = genes[pos:]
        karyo[cid] = genes[:pos]
        if centros:
            centros[cid] = pos
            centros[new_id] = 0
        self.events.append({"node": node_name, "type": "fission",
                            "chroms": [cid, new_id], "pos": "split {} at {}".format(cid, pos)})
        return True

    def _apply_gene_gain(self, karyo, node_name, centros=None):
        if not karyo:
            return False
        cid = self.rng.choice(list(karyo.keys()))
        genes = karyo[cid]
        new_gid = self.tracker.new_gain_gene(species=None)
        orient = self.rng.choice(["+", "-"])
        pos = self.rng.randint(0, len(genes))
        genes.insert(pos, (new_gid, orient))
        if centros and cid in centros:
            c = centros[cid]
            if c > pos:
                centros[cid] = c + 1
        self.events.append({"node": node_name, "type": "gene_gain",
                            "chrom": cid, "gene": new_gid, "pos": pos})
        return True

    def _apply_gene_loss(self, karyo, node_name, centros=None):
        if not karyo:
            return False
        cids = [c for c in karyo if len(karyo[c]) > 0]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        idx = self.rng.randint(0, len(genes) - 1)
        lost_gid, lost_orient = genes.pop(idx)
        self.tracker.gene_ancestor.pop(lost_gid, None)
        self.tracker.gene_species.pop(lost_gid, None)
        self.tracker.gene_wgd_copies.pop(lost_gid, None)
        self.tracker.gene_is_gain.pop(lost_gid, None)
        self.tracker.gene_dup_of.pop(lost_gid, None)
        self.tracker.gene_dup_type.pop(lost_gid, None)
        if centros and cid in centros:
            c = centros[cid]
            if c > idx:
                centros[cid] = c - 1
        if not genes:
            del karyo[cid]
            if centros:
                centros.pop(cid, None)
        self.events.append({"node": node_name, "type": "gene_loss",
                            "chrom": cid, "gene": lost_gid, "pos": idx})
        return True

    def _sample_seg_len(self, max_len):
        """Geometric-like distribution: shorter segments more likely."""
        seg_len = 1
        while seg_len < max_len and self.rng.random() < 0.5:
            seg_len += 1
        return seg_len

    def _apply_tandem_dup(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 1]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        max_len = min(10, len(genes))
        seg_len = self._sample_seg_len(max_len)
        start = self.rng.randint(0, len(genes) - seg_len)
        segment = genes[start:start + seg_len]
        dup_segment = []
        for gid, orient in segment:
            new_gid = self.tracker.new_duplicate_gene(gid, "tandem")
            dup_segment.append((new_gid, orient))
        insert_pos = start + seg_len
        genes[insert_pos:insert_pos] = dup_segment
        if centros and cid in centros:
            c = centros[cid]
            if c > insert_pos:
                centros[cid] = c + seg_len
        self.events.append({"node": node_name, "type": "tandem_dup",
                            "chrom": cid, "pos": start,
                            "len": seg_len, "genes": [g for g, _ in segment]})
        return True

    def _apply_dispersed_dup(self, karyo, node_name, centros=None):
        if not karyo:
            return False
        src_cid = self.rng.choice(list(karyo.keys()))
        src_genes = karyo[src_cid]
        if not src_genes:
            return False
        max_len = min(10, len(src_genes))
        seg_len = self._sample_seg_len(max_len)
        start = self.rng.randint(0, len(src_genes) - seg_len)
        segment = src_genes[start:start + seg_len]
        dup_segment = []
        for gid, orient in segment:
            new_gid = self.tracker.new_duplicate_gene(gid, "dispersed")
            dup_segment.append((new_gid, orient))
        tgt_cid = self.rng.choice(list(karyo.keys()))
        tgt_genes = karyo[tgt_cid]
        tgt_pos = self.rng.randint(0, len(tgt_genes))
        tgt_genes[tgt_pos:tgt_pos] = dup_segment
        if centros and tgt_cid in centros:
            c = centros[tgt_cid]
            if c > tgt_pos:
                centros[tgt_cid] = c + seg_len
        self.events.append({"node": node_name, "type": "dispersed_dup",
                            "source_chrom": src_cid, "target_chrom": tgt_cid,
                            "target_pos": tgt_pos, "len": seg_len,
                            "genes": [g for g, _ in segment]})
        return True

    def _apply_unidir_trans(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 3]
        if len(cids) < 2:
            return False
        src_cid, tgt_cid = self.rng.sample(cids, 2)
        src_genes = karyo[src_cid]
        max_len = min(10, len(src_genes) - 2)
        seg_len = self._sample_seg_len(max_len)
        start = self.rng.randint(1, len(src_genes) - seg_len - 1)
        if centros and src_cid in centros:
            c = centros[src_cid]
            if start <= c < start + seg_len:
                return False
        segment = src_genes[start:start + seg_len]
        if self.rng.random() < 0.5:
            segment = reverse_segment(segment)
        del src_genes[start:start + seg_len]
        if centros and src_cid in centros:
            c = centros[src_cid]
            if c > start + seg_len:
                centros[src_cid] = c - seg_len
        tgt_genes = karyo[tgt_cid]
        tgt_pos = self.rng.randint(0, len(tgt_genes))
        tgt_genes[tgt_pos:tgt_pos] = segment
        if centros and tgt_cid in centros:
            c = centros[tgt_cid]
            if c > tgt_pos:
                centros[tgt_cid] = c + seg_len
        self.events.append({"node": node_name, "type": "unidir_trans",
                            "source_chrom": src_cid, "target_chrom": tgt_cid,
                            "seg_len": seg_len, "pos": start})
        return True

    def _apply_seg_deletion(self, karyo, node_name, centros=None):
        cids = [c for c in karyo if len(karyo[c]) >= 5]
        if not cids:
            return False
        cid = self.rng.choice(cids)
        genes = karyo[cid]
        max_len = min(50, len(genes) // 20)
        if max_len < 2:
            return False
        seg_len = self.rng.randint(2, max_len)
        start = self.rng.randint(0, len(genes) - seg_len)
        del genes[start:start + seg_len]
        if centros and cid in centros:
            c = centros[cid]
            if c > start + seg_len:
                centros[cid] = c - seg_len
            elif c > start:
                centros[cid] = start
        if not genes:
            del karyo[cid]
            if centros:
                centros.pop(cid, None)
        self.events.append({"node": node_name, "type": "seg_deletion",
                            "chrom": cid, "pos": "{}-{}".format(start, start + seg_len),
                            "len": seg_len})
        return True

    def _apply_seg_duplication(self, karyo, node_name, centros=None):
        if not karyo:
            return False
        src_cid = self.rng.choice(list(karyo.keys()))
        src_genes = karyo[src_cid]
        if len(src_genes) < 20:
            return False
        max_len = min(30, len(src_genes) // 10)
        if max_len < 10:
            return False
        seg_len = self.rng.randint(10, max_len)
        start = self.rng.randint(0, len(src_genes) - seg_len)
        segment = src_genes[start:start + seg_len]
        if self.rng.choice([True, False]):
            tgt_cid = src_cid
            tgt_genes = src_genes
            tgt_pos = start + seg_len
        else:
            tgt_cid = self.rng.choice(list(karyo.keys()))
            tgt_genes = karyo[tgt_cid]
            tgt_pos = self.rng.randint(0, len(tgt_genes))
        dup_segment = []
        for gid, orient in segment:
            new_gid = self.tracker.new_duplicate_gene(gid, "seg_dup")
            dup_segment.append((new_gid, orient))
        tgt_genes[tgt_pos:tgt_pos] = dup_segment
        if centros and tgt_cid in centros:
            c = centros[tgt_cid]
            if c > tgt_pos:
                centros[tgt_cid] = c + seg_len
        self.events.append({"node": node_name, "type": "seg_duplication",
                            "source_chrom": src_cid, "target_chrom": tgt_cid,
                            "pos": start, "len": seg_len})
        return True

    def _apply_chromothripsis(self, karyo, node_name, centros=None):
        if not karyo:
            return False
        cid = self.rng.choice(list(karyo.keys()))
        genes = karyo[cid]
        if len(genes) < 10:
            return False
        n_breaks = self.rng.randint(3, min(10, len(genes) - 1))
        break_points = sorted(self.rng.sample(range(1, len(genes)), n_breaks))
        segments = []
        prev = 0
        for bp in break_points:
            segments.append(genes[prev:bp])
            prev = bp
        segments.append(genes[prev:])
        kept = [seg for seg in segments if self.rng.random() > 0.3]
        if not kept:
            return False
        self.rng.shuffle(kept)
        new_genes = []
        for seg in kept:
            if self.rng.choice([True, False]):
                seg = reverse_segment(seg)
            new_genes.extend(seg)
        if not new_genes:
            return False
        karyo[cid] = new_genes
        if centros and cid in centros:
            centros[cid] = self.rng.randint(2, len(new_genes) - 2)
        self.events.append({"node": node_name, "type": "chromothripsis",
                            "chrom": cid, "breaks": n_breaks,
                            "kept": len(kept), "lost": len(segments) - len(kept)})
        return True

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def _print_summary(self):
        print("\n=== Simulation Summary ===")
        print("Ancestral genes: {}".format(len(self.tracker.ancestor_genes)))
        print("Total genes created: {}".format(self.tracker.total_genes()))
        print("Gene gains: {}".format(len(self.tracker.gain_genes)))
        print("Species simulated: {}".format(len(self.leaf_karyotypes)))
        for sp in sorted(self.leaf_karyotypes):
            karyo = self.leaf_karyotypes[sp]
            n_chr = len(karyo)
            n_genes = sum(len(g) for g in karyo.values())
            p = self.wgd_map.get(sp, 1)
            p_str = " [p={}]".format(p) if p > 1 else ""
            print("  {}: {} chromosomes, {} genes{}".format(sp, n_chr, n_genes, p_str))
        event_counts = defaultdict(int)
        gene_totals = defaultdict(int)
        for e in self.events:
            etype = e["type"]
            event_counts[etype] += 1
            if "len" in e:
                gene_totals[etype] += e["len"]
            elif "seg_len" in e:
                gene_totals[etype] += e["seg_len"]
            elif etype == "gene_gain" or etype == "gene_loss":
                gene_totals[etype] += 1
            elif etype == "fractionation":
                gene_totals[etype] += e.get("lost", 0)
            elif etype == "chromothripsis":
                gene_totals[etype] += e.get("lost", 0)
            elif etype == "inversion":
                try:
                    p1, p2 = e.get("pos", "0-0").split("-")
                    gene_totals[etype] += int(p2) - int(p1)
                except (ValueError, IndexError):
                    pass
        print("Events:")
        order = [
            "WGD", "fission", "EEJ", "NCF",
            "RT", "inversion", "chromothripsis",
            "seg_deletion", "seg_duplication",
            "tandem_dup", "dispersed_dup", "unidir_trans",
			"gene_gain", "gene_loss", "fractionation",
        ]
        printed = set()
        for etype in order:
            if etype in event_counts:
                gc = gene_totals.get(etype, 0)
                if gc:
                    print("  {}: {} ({} genes)".format(etype, event_counts[etype], gc))
                else:
                    print("  {}: {}".format(etype, event_counts[etype]))
                printed.add(etype)
        for etype, count in sorted(event_counts.items()):
            if etype not in printed and etype != "rearrangements":
                gc = gene_totals.get(etype, 0)
                if gc:
                    print("  {}: {} ({} genes)".format(etype, count, gc))
                else:
                    print("  {}: {}".format(etype, count))

    # ==================================================================
    # Output generation
    # ==================================================================
    def generate_outputs(self, outdir, tree, ploidy_map, orig_nw):
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        species = sorted(self.leaf_karyotypes.keys())
        print("\nGenerating output files in {} ...".format(outdir))
        output_files = []

        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    self.tracker.set_species(gid, sp)

        # Write tree first
        self._write_tree(outdir, orig_nw, tree)
        output_files.append(os.path.join(outdir, "species_tree.nwk"))

        self._write_leaf_gff(outdir, species)
        output_files.append(os.path.join(outdir, "all_species_gene.gff"))
        self._write_ancestor_gff(outdir, tree)
        output_files.append(os.path.join(outdir, "ancestors_gene.gff"))
        self._write_og(outdir, species)
        output_files.append(os.path.join(outdir, "ortholog_groups.txt"))
        self._write_orthologs(outdir, species)
        output_files.append(os.path.join(outdir, "ortholog_pairs.txt"))
        self._write_events(outdir, tree)
        output_files.append(os.path.join(outdir, "events.tsv"))
        self._write_true_hogs(outdir, tree, species)
        output_files.append(os.path.join(outdir, "true_hogs.tsv"))
        self._write_ancestor_karyotypes(outdir, tree)
        output_files.append(os.path.join(outdir, "ancestors_karyotypes.txt"))
        self._write_gene_ancestor_map(outdir, tree)
        output_files.append(os.path.join(outdir, "gene_ancestor_map.tsv"))
        draw_ancestor_dotplots(
            tree, self.all_node_karyotypes, outdir,
            ancestor_fn=self.tracker.get_ancestor)
        # Per-branch dotplots (parent → child)
        branch_dir = os.path.join(outdir, "branch_dotplots")
        draw_branch_dotplots(
            tree, self.all_node_karyotypes, branch_dir,
            gene_parent=self.tracker.gene_parent)
        # Sister dotplots (sibling ↔ sibling)
        sister_dir = os.path.join(outdir, "sister_dotplots")
        draw_sister_dotplots(
            tree, self.all_node_karyotypes, sister_dir)
        # Collect dotplot PNGs
        root_name = None
        for node in tree.traverse("preorder"):
            if node.is_root():
                root_name = node.name
                break
        for node in tree.traverse("preorder"):
            name = node.name
            if name == root_name or name not in self.all_node_karyotypes:
                continue
            output_files.append(os.path.join(outdir, "dotplot_{}.png".format(name)))
        # Collect branch dotplot PNGs
        for node in tree.traverse("preorder"):
            if node.is_root() or node.up is None:
                continue
            parent_key = node.up.name
            child_key = node.name
            if parent_key in self.all_node_karyotypes and child_key in self.all_node_karyotypes:
                parent_label = parent_key if parent_key else "root"
                child_label = child_key if child_key else "unnamed"
                output_files.append(os.path.join(
                    branch_dir, "dotplot_{}_{}.png".format(parent_label, child_label)))
        # Collect sister dotplot PNGs
        for node in tree.traverse("preorder"):
            if node.is_leaf():
                continue
            children = [c for c in node.children if c.name in self.all_node_karyotypes]
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    c1, c2 = children[i].name, children[j].name
                    output_files.append(os.path.join(
                        sister_dir, "dotplot_{}_vs_{}.png".format(c1, c2)))

        print("Done! Output files:")
        for f in sorted(output_files):
            print("  {}".format(f))

    # ------------------------------------------------------------------
    # Tree output — with [p=N] for simulated WGD
    # ------------------------------------------------------------------
    def _write_tree(self, outdir, orig_nw, tree):
        """Write tree file with [p=N] annotations for WGD species."""
        path = os.path.join(outdir, "species_tree.nwk")

        # Build the tree Newick with [p=N] annotations
        # Use ete3 to write, then add [p=N] to leaf names
        nw = tree.write(format=1).strip()
        if not nw.endswith(';'):
            nw += ';'

        # Add [p=N] annotations to species that had WGD
        for sp, factor in self.wgd_map.items():
            if factor > 1:
                # Replace SpeciesName: with SpeciesName[p=2]: etc.
                nw = re.sub(
                    r'(\b{})(:|,)'.format(re.escape(sp)),
                    r'\1[p={}]\2'.format(factor),
                    nw)

        with open(path, 'w') as f:
            f.write(nw + '\n')

    # ------------------------------------------------------------------
    # Leaf GFF
    # ------------------------------------------------------------------
    def _write_leaf_gff(self, outdir, species):
        path = os.path.join(outdir, "all_species_gene.gff")
        with open(path, 'w') as f:
            for sp in species:
                karyo = self.leaf_karyotypes[sp]
                sorted_cids = self._sorted_chr_ids(karyo)
                for cid in sorted_cids:
                    chrom_id = "{}|{}".format(sp, cid)
                    for i, (gid, orient) in enumerate(karyo[cid]):
                        start = (i + 1) * 100 - 99
                        end = (i + 1) * 100
                        f.write("{}\t{}|{}\t{}\t{}\t{}\n".format(
                            chrom_id, sp, gid, start, end, orient))

    def _write_ancestor_gff(self, outdir, tree):
        path = os.path.join(outdir, "ancestors_gene.gff")
        with open(path, 'w') as f:
            for node in tree.traverse("preorder"):
                if node.name not in self.all_node_karyotypes:
                    continue
                karyo = self.all_node_karyotypes[node.name]
                sorted_cids = self._sorted_chr_ids(karyo)
                for cid in sorted_cids:
                    chrom_id = "{}|{}".format(node.name, cid)
                    for i, (gid, orient) in enumerate(karyo[cid]):
                        start = (i + 1) * 100 - 99
                        end = (i + 1) * 100
                        f.write("{}\t{}|{}\t{}\t{}\t{}\n".format(
                            chrom_id, node.name, gid, start, end, orient))

    def _write_ancestor_karyotypes(self, outdir, tree):
        path = os.path.join(outdir, "ancestors_karyotypes.txt")
        with open(path, 'w') as f:
            for node in tree.traverse("preorder"):
                if node.name not in self.all_node_karyotypes:
                    continue
                karyo = self.all_node_karyotypes[node.name]
                f.write(">{}\t{} chromosomes\n".format(node.name, len(karyo)))
                sorted_cids = self._sorted_chr_ids(karyo)
                for cid in sorted_cids:
                    chrom_id = "{}|{}".format(node.name, cid)
                    genes = karyo[cid]
                    gene_str = " ".join(
                        ("-" + g if o == "-" else g) for g, o in genes)
                    f.write("{}\t{}\n".format(chrom_id, gene_str))
            # Write pre-WGD snapshots (not in tree, but in all_node_karyotypes)
            for name, karyo in sorted(self.all_node_karyotypes.items()):
                if not name.endswith("_preWGD"):
                    continue
                f.write(">{}\t{} chromosomes\n".format(name, len(karyo)))
                sorted_cids = self._sorted_chr_ids(karyo)
                for cid in sorted_cids:
                    chrom_id = "{}|{}".format(name, cid)
                    genes = karyo[cid]
                    gene_str = " ".join(
                        ("-" + g if o == "-" else g) for g, o in genes)
                    f.write("{}\t{}\n".format(chrom_id, gene_str))

    def _write_gene_ancestor_map(self, outdir, tree):
        """Write gene → ancestor (HOG) mapping.

        Output: gene_ancestor_map.tsv
        Columns: gene, species, ancestor, wgd_copies, chromosome

        The 'ancestor' column is the true HOG identifier — genes sharing
        the same ancestor belong to the same ortholog group (SOG).
        """
        path = os.path.join(outdir, "gene_ancestor_map.tsv")
        with open(path, 'w') as f:
            f.write("gene\tspecies\tancestor\twgd_copies\tchromosome\n")
            for node in tree.traverse("preorder"):
                name = node.name
                if name not in self.all_node_karyotypes:
                    continue
                karyo = self.all_node_karyotypes[name]
                for cid, genes in karyo.items():
                    for gid, orient in genes:
                        anc = self.tracker.get_ancestor(gid)
                        wgd = self.tracker.gene_wgd_copies.get(gid, {})
                        wgd_str = ";".join(
                            "{}:{}".format(k, v) for k, v in sorted(wgd.items())
                        ) if wgd else ""
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(
                            gid, name, anc, wgd_str, cid))

    def _write_og(self, outdir, species):
        anc_to_genes = defaultdict(list)
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    anc = self.tracker.get_ancestor(gid)
                    full_name = "{}|{}".format(sp, gid)
                    anc_to_genes[anc].append(full_name)
        path = os.path.join(outdir, "ortholog_groups.txt")
        with open(path, 'w') as f:
            for i, (anc, gene_sp_list) in enumerate(sorted(anc_to_genes.items()), 1):
                gene_names = sorted(gene_sp_list)
                f.write("SOG{}: {}\n".format(i, " ".join(gene_names)))

    def _write_orthologs(self, outdir, species):
        sp_gene_groups = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_gene_groups[sp][(anc, wgd_key)].append(full_name)
        path = os.path.join(outdir, "ortholog_pairs.txt")
        with open(path, 'w') as f:
            sp_list = sorted(species)
            for i in range(len(sp_list)):
                for j in range(i + 1, len(sp_list)):
                    sp1, sp2 = sp_list[i], sp_list[j]
                    groups1 = sp_gene_groups[sp1]
                    groups2 = sp_gene_groups[sp2]
                    common_keys = set(groups1.keys()) & set(groups2.keys())
                    for key in common_keys:
                        for g1 in sorted(groups1[key]):
                            for g2 in sorted(groups2[key]):
                                f.write("{}\t{}\n".format(g1, g2))
            self._write_cross_wgd_orthologs(f, species, sp_gene_groups)

    def _write_cross_wgd_orthologs(self, fout, species, sp_gene_groups):
        sp_anc_genes = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_genes[sp][anc].append((full_name, wgd_key))
        sp_list = sorted(species)
        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                common_ancs = set(sp_anc_genes[sp1].keys()) & set(sp_anc_genes[sp2].keys())
                for anc in common_ancs:
                    genes1 = sp_anc_genes[sp1][anc]
                    genes2 = sp_anc_genes[sp2][anc]
                    for g1, wk1 in genes1:
                        for g2, wk2 in genes2:
                            if wk1 == wk2:
                                continue
                            shared_wgd_nodes = set(k for k, _ in wk1) & set(k for k, _ in wk2)
                            compatible = all(
                                dict(wk1).get(n) == dict(wk2).get(n)
                                for n in shared_wgd_nodes)
                            if compatible:
                                fout.write("{}\t{}\n".format(min(g1, g2), max(g1, g2)))

    def _write_events(self, outdir, tree, path="events.tsv"):
        """Write events in unified branch-level format.

        Format: branch, event_type, genes, ancestors, chroms, desc, support
        where branch = parent->child identifier,
              ancestors = comma-separated ancestor (HOG) IDs for each gene.

        For WGD nodes, an extra row is emitted for the parent→preWGD branch
        so that the reconstruction's pre-WGD node has a corresponding truth.
        """
        # Build parent_of mapping from tree
        parent_of = {}
        for node in tree.traverse("preorder"):
            if not node.is_root():
                parent_of[node.name] = node.up.name

        # Collect WGD nodes for pre-WGD branch emission
        wgd_nodes = set()
        for e in self.events:
            if e['type'] == 'WGD':
                wgd_nodes.add(e['node'])

        evt_path = os.path.join(outdir, path)
        with open(evt_path, 'w') as f:
            f.write("branch\tevent_type\tgenes\tancestors\tchroms\tdesc\tsupport\n")
            # Emit pre-WGD branch placeholder for each WGD node
            for wgd_node in sorted(wgd_nodes):
                parent = parent_of.get(wgd_node, '?')
                pre_wgd = "{}_preWGD".format(wgd_node)
                f.write("{}-{}\t{}\t\t\t\tpre-WGD snapshot\n".format(
                    parent, pre_wgd, "pre_wgd_reference"))

            for e in self.events:
                if e['type'] in ('rearrangements',):
                    continue  # skip summary events
                # branch = parent - child
                child = e['node']
                parent = parent_of.get(child, '?')
                branch = "{}-{}".format(parent, child)
                # canonicalize event type
                from .takr_events import canonicalize_event_type
                etype = canonicalize_event_type(e['type'])
                # extract fields
                genes = self._fmt_genes(e)
                # build ancestor (HOG) string for involved genes
                ancestors = self._fmt_ancestors(e)
                chroms = e.get('chroms', e.get('chrom', ''))
                if isinstance(chroms, list):
                    chroms = ','.join(str(x) for x in chroms)
                desc = e.get('desc', e.get('pos', ''))
                if not desc:
                    desc = "{}: {}".format(etype, chroms)
                support = e.get('support', 1)
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    branch, etype, genes, ancestors, chroms, desc, support))

    def _fmt_genes(self, e):
        """Extract genes string from event dict."""
        for key in ('genes', 'gene'):
            vals = e.get(key)
            if vals:
                if isinstance(vals, list):
                    return ','.join(str(v) for v in vals)
                return str(vals)
        return ''

    def _fmt_ancestors(self, e):
        """Map genes in event to their ancestor (HOG) IDs."""
        genes = e.get('genes') or e.get('gene')
        if not genes:
            # fractionation events have 'gene' (singular)
            gid = e.get('gene')
            if gid:
                anc = self.tracker.get_ancestor(gid)
                return str(anc) if anc else ''
            return ''
        if isinstance(genes, list):
            ancs = []
            for g in genes:
                anc = self.tracker.get_ancestor(g)
                ancs.append(str(anc) if anc else '')
            return ','.join(ancs)
        # single gene string
        anc = self.tracker.get_ancestor(str(genes))
        return str(anc) if anc else ''

    def _write_true_hogs(self, outdir, tree, species):
        tree_copy = tree.copy("deepcopy")
        number_tree_nodes(tree_copy)

        G = nx.Graph()
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    full_name = "{}|{}".format(sp, gid)
                    G.add_node(full_name)

        sp_anc_genes = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_genes[sp][(anc, wgd_key)].append(full_name)

        sp_list = sorted(species)
        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                groups1 = sp_anc_genes[sp1]
                groups2 = sp_anc_genes[sp2]
                common_keys = set(groups1.keys()) & set(groups2.keys())
                for key in common_keys:
                    for g1 in groups1[key]:
                        for g2 in groups2[key]:
                            G.add_edge(g1, g2)

        sp_anc_only = defaultdict(lambda: defaultdict(list))
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    if self.tracker.gene_is_gain.get(gid, False):
                        continue
                    anc = self.tracker.get_ancestor(gid)
                    wgd = self.tracker.gene_wgd_copies.get(gid, {})
                    wgd_key = tuple(sorted(wgd.items()))
                    full_name = "{}|{}".format(sp, gid)
                    sp_anc_only[sp][anc].append((full_name, wgd_key))

        for i in range(len(sp_list)):
            for j in range(i + 1, len(sp_list)):
                sp1, sp2 = sp_list[i], sp_list[j]
                common_ancs = set(sp_anc_only[sp1].keys()) & set(sp_anc_only[sp2].keys())
                for anc in common_ancs:
                    genes1 = sp_anc_only[sp1][anc]
                    genes2 = sp_anc_only[sp2][anc]
                    for g1, wk1 in genes1:
                        for g2, wk2 in genes2:
                            if wk1 == wk2:
                                continue
                            shared = set(k for k, _ in wk1) & set(k for k, _ in wk2)
                            if all(dict(wk1).get(n) == dict(wk2).get(n) for n in shared):
                                G.add_edge(g1, g2)

        og_genes = defaultdict(list)
        for sp in species:
            for cid, genes in self.leaf_karyotypes[sp].items():
                for gid, orient in genes:
                    anc = self.tracker.get_ancestor(gid)
                    full_name = "{}|{}".format(sp, gid)
                    og_genes[anc].append(full_name)
        for gid in self.tracker.gain_genes:
            sp = self.tracker.gene_species.get(gid, "unknown")
            full_name = "{}|{}".format(sp, gid)
            og_genes[gid].append(full_name)

        all_hogs = []
        node_species = {}
        for node in tree_copy.traverse("postorder"):
            if node.is_leaf():
                node_species[node.name] = {node.name}
            else:
                node_species[node.name] = set()
                for child in node.children:
                    node_species[node.name] |= node_species[child.name]

        og_counter = 0
        for anc in sorted(og_genes.keys()):
            og_counter += 1
            og_id = "SOG{}".format(og_counter)
            og_gene_set = set(og_genes[anc])

            node_hogs = {}
            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                sp_set = node_species[nid]
                node_genes = set()
                for g in og_gene_set:
                    sp = g.split("|")[0]
                    if sp in sp_set:
                        node_genes.add(g)
                if not node_genes:
                    continue
                subG = G.subgraph(node_genes)
                ccs = list(nx.connected_components(subG))
                hog_list = []
                for idx, cc in enumerate(sorted(ccs, key=lambda x: sorted(x)[0] if x else "")):
                    hog_id = "{}.{}.hog{}".format(og_id, nid, idx)
                    hog_list.append((hog_id, cc))
                node_hogs[nid] = hog_list

            hog_parent = {}
            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                if nid not in node_hogs:
                    continue
                if node.up is None:
                    for hog_id, cc in node_hogs[nid]:
                        hog_parent[hog_id] = "Root"
                else:
                    parent_nid = node.up.name
                    if parent_nid not in node_hogs:
                        for hog_id, cc in node_hogs[nid]:
                            hog_parent[hog_id] = "Root"
                    else:
                        for hog_id, cc in node_hogs[nid]:
                            found_parent = False
                            for parent_hog_id, parent_cc in node_hogs[parent_nid]:
                                if cc.issubset(parent_cc):
                                    hog_parent[hog_id] = parent_hog_id
                                    found_parent = True
                                    break
                            if not found_parent:
                                hog_parent[hog_id] = "Root"

            for node in tree_copy.traverse("postorder"):
                if node.is_leaf():
                    continue
                nid = node.name
                if nid not in node_hogs:
                    continue
                for hog_id, cc in node_hogs[nid]:
                    all_hogs.append({
                        "hog_id": hog_id, "og_id": og_id,
                        "node_id": nid, "parent": hog_parent.get(hog_id, "Root"),
                        "genes": sorted(cc)})

        path = os.path.join(outdir, "true_hogs.tsv")
        with open(path, 'w') as f:
            f.write("{}\n".format("\t".join(["HOG", "OG", "Tree_Node", "Parent", "Genes"])))
            for hog in all_hogs:
                gene_str = " ".join(hog["genes"])
                f.write("{}\n".format("\t".join([
                    hog["hog_id"], hog["og_id"], hog["node_id"],
                    hog["parent"], gene_str])))

    # ------------------------------------------------------------------
    # Dotplot drawing
    # ------------------------------------------------------------------
    # Helper: sort chromosome IDs
    # ------------------------------------------------------------------
    def _sorted_chr_ids(self, karyo):
        def sort_key(c):
            m = re.search(r'(\d+)', c)
            return (0, int(m.group(1))) if m else (1, c)
        return sorted(karyo.keys(), key=sort_key)


# ============================================================================
# CLI argument definition
# ============================================================================

def sim_args(parser):
    """Define arguments for the sim subcommand."""
    g_input = parser.add_argument_group('Input')
    g_input.add_argument('-t', '-tree', type=str, default=None,
                         dest='tree', metavar='FILE',
                         help='Newick tree with branch lengths')
    g_input.add_argument('-n', '--num-species', type=int, default=None,
                         dest='num_species', metavar='INT',
                         help='Number of species (auto-generate tree)')

    parser.add_argument('-o', '-outdir', type=str, default='./',
                        dest='outdir', metavar='DIR',
                        help='Output directory [default=%(default)s]')

    g_anc = parser.add_argument_group('Ancestor genome')
    g_anc.add_argument('--num-chroms', type=int, default=8, metavar='INT',
                       help='Number of ancestral chromosomes [default=%(default)s]')
    g_anc.add_argument('--min-genes', type=int, default=200, metavar='INT',
                       help='Min genes per ancestral chromosome [default=%(default)s]')
    g_anc.add_argument('--max-genes', type=int, default=1000, metavar='INT',
                       help='Max genes per ancestral chromosome [default=%(default)s]')

    g_rate = parser.add_argument_group('Rates (per unit branch length)')
    g_rate.add_argument('--inv-rate', type=float, metavar='FLOAT', default=10.0,
                        help='Inversion rate [default=%(default)s]')
    g_rate.add_argument('--rt-rate', type=float, metavar='FLOAT', default=4.0,
                        help='Reciprocal translocation rate [default=%(default)s]')
    g_rate.add_argument('--ncf-rate', type=float, metavar='FLOAT', default=2.0,
                        help='Non-centric fusion rate [default=%(default)s]')
    g_rate.add_argument('--eej-rate', type=float, metavar='FLOAT', default=2.0,
                        help='End-end joining rate [default=%(default)s]')
    g_rate.add_argument('--fission-rate', type=float, metavar='FLOAT', default=0.02,
                        help='Fission rate (very low) [default=%(default)s]')
    g_rate.add_argument('--wgd-rate', type=float, metavar='FLOAT', default=0.2,
                        help='WGD rate (2x most common) [default=%(default)s]')
    g_rate.add_argument('--gene-gain-rate', type=float, metavar='FLOAT', default=0.5,
                        help='Gene gain rate [default=%(default)s]')
    g_rate.add_argument('--gene-loss-rate', type=float, metavar='FLOAT', default=0.5,
                        help='Gene loss rate [default=%(default)s]')
    g_rate.add_argument('--tandem-dup-rate', type=float, metavar='FLOAT', default=10.0,
                        help='Tandem duplication rate [default=%(default)s]')
    g_rate.add_argument('--dispersed-dup-rate', type=float, metavar='FLOAT', default=2.0,
                        help='Dispersed duplication rate [default=%(default)s]')
    g_rate.add_argument('--unidir-trans-rate', type=float, metavar='FLOAT', default=1.0,
                        help='Unidirectional translocation rate [default=%(default)s]')
    g_rate.add_argument('--frac-rate', type=float, metavar='FLOAT', default=1000,
                        help='Fractionation rate after WGD (Poisson rate per gene) [default=%(default)s]')
    g_rate.add_argument('--seg-del-rate', type=float, metavar='FLOAT', default=0.1,
                        help='Segmental deletion rate [default=%(default)s]')
    g_rate.add_argument('--seg-dup-rate', type=float, metavar='FLOAT', default=0.2,
                        help='Segmental duplication rate [default=%(default)s]')
    g_rate.add_argument('--chromothripsis-rate', type=float, metavar='FLOAT', default=0.01,
                        help='Chromothripsis rate [default=%(default)s]')
    g_rate.add_argument('--scale', type=float, metavar='FLOAT', default=1.0,
                        help='Global rate scaling factor [default=%(default)s]')

    parser.add_argument('--seed', type=int, default=42, metavar='INT',
                        help='Random seed [default=%(default)s]')


def xmain(**kargs):
    """Main entry point for the sim subcommand."""
    seed = kargs.get('seed', 42)
    rng = random.Random(seed)

    tree_file = kargs.get('tree', None)
    num_species = kargs.get('num_species', None)

    if tree_file:
        wgd_rate = kargs.get('wgd_rate', 0.4)
        tree, ploidy_map, orig_nw = parse_tree(tree_file, wgd_rate)
        if ploidy_map:
            print("Ploidy annotations found: {}".format(ploidy_map))
    elif num_species:
        if num_species < 2:
            print("ERROR: --num-species must be >= 2", file=sys.stderr)
            sys.exit(1)
        tree = generate_tree(num_species, rng)
        # Renumber internal nodes
        existing = set(n.name for n in tree.traverse() if n.name)
        counter = 0
        for node in tree.traverse("preorder"):
            if not node.name or node.name.strip() == "":
                while True:
                    candidate = "N{}".format(counter)
                    counter += 1
                    if candidate not in existing:
                        node.name = candidate
                        existing.add(candidate)
                        break
        ploidy_map = {}
        orig_nw = tree.write(format=1)
        print("Generated tree: {}".format(orig_nw))
    else:
        print("ERROR: either -t or --num-species is required", file=sys.stderr)
        sys.exit(1)

    if kargs.get('wgd_rate', 0.4) > 0:
        print("WGD rate: {} (rate-based, factor sampled: 2x most common)".format(
            kargs['wgd_rate']))

    sim = EvolutionSimulator(
        seed=seed,
        num_chroms=kargs.get('num_chroms', 8),
        min_genes=kargs.get('min_genes', 200),
        max_genes=kargs.get('max_genes', 1000),
        inv_rate=kargs.get('inv_rate', 10.0),
        rt_rate=kargs.get('rt_rate', 3.0),
        ncf_rate=kargs.get('ncf_rate', 1.5),
        eej_rate=kargs.get('eej_rate', 2.0),
        fission_rate=kargs.get('fission_rate', 0.02),
        wgd_rate=kargs.get('wgd_rate', 0.4),
        gene_gain_rate=kargs.get('gene_gain_rate', 0.5),
        tandem_dup_rate=kargs.get('tandem_dup_rate', 2.0),
        dispersed_dup_rate=kargs.get('dispersed_dup_rate', 0.3),
        unidir_trans_rate=kargs.get('unidir_trans_rate', 1.0),
        frac_rate=kargs.get('frac_rate', 5.0),
        seg_del_rate=kargs.get('seg_del_rate', 0.5),
        seg_dup_rate=kargs.get('seg_dup_rate', 0.3),
        chromothripsis_rate=kargs.get('chromothripsis_rate', 0.05),
        scale=kargs.get('scale', 1.0),
    )

    sim.run(tree, ploidy_map)
    sim.generate_outputs(kargs['outdir'], tree, ploidy_map, orig_nw)


def main():
    parser = argparse.ArgumentParser(
        description="Chromosome Rearrangement Evolution Simulator (AK.py format)",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sim_args(parser)
    args = parser.parse_args()
    if not args.tree and not args.num_species:
        print("ERROR: either -t or --num-species is required", file=sys.stderr)
        sys.exit(1)
    xmain(**vars(args))


if __name__ == "__main__":
    main()
