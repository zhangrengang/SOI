#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dotplot drawing utilities for ancestral karyotype reconstruction.

Provides reusable functions to draw synteny dotplots between
an ancestral genome and its descendants, or between any two karyotypes.
"""

import os
import re

from . import style  # global mpl settings


def _sorted_chr_ids(karyo):
    """Sort chromosome IDs numerically if they contain digits."""
    def sort_key(c):
        m = re.search(r'(\d+)', c)
        return (0, int(m.group(1))) if m else (1, c)
    return sorted(karyo.keys(), key=sort_key)


def _gene_id(g):
    """Extract gene ID string from a tuple (gid, orient) or plain string."""
    if isinstance(g, tuple):
        return g[0]
    return g


def _chr_boundaries(ax, chr_list, karyo, axis='x'):
    """Draw chromosome boundaries and return label positions."""
    pos = 0
    labels = []
    for cid in chr_list:
        n_genes = len(karyo[cid])
        if pos > 0:
            if axis == 'x':
                ax.axvline(pos, color='black', linewidth=1.5, alpha=0.5)
            else:
                ax.axhline(pos, color='black', linewidth=1.5, alpha=0.5)
        labels.append((pos + n_genes / 2.0, cid))
        pos += n_genes
    return labels, pos


def draw_dotplot(ref_karyo, query_karyo, outpath,
                 ref_name="Reference", query_name="Query",
                 gene_mapping=None, ref_chrom_of=None, dpi=200):
    """Draw a dotplot comparing query karyotype against reference karyotype.

    Parameters
    ----------
    ref_karyo : dict
        {chrom_id: [gene_id, ...]} reference karyotype.
    query_karyo : dict
        {chrom_id: [gene_id, ...]} query karyotype.
    outpath : str
        Output PNG file path.
    ref_name, query_name : str
        Labels for the axes.
    gene_mapping : callable or dict, optional
        Maps a query gene to a reference gene. If None, identity mapping
        is used.
    ref_chrom_of : dict, optional
        Maps reference gene_id to its chromosome id for coloring.
        If None, colors are derived from ref_karyo.
    dpi : int
        Output resolution.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        print("WARNING: matplotlib not available, skipping dotplot generation")
        return

    if gene_mapping is None:
        gene_mapping = lambda g: g
    elif isinstance(gene_mapping, dict):
        _map = gene_mapping
        gene_mapping = lambda g: _map.get(g, g)

    # Build reference gene order and chromosome colors
    ref_order = {}
    pos = 0
    sorted_ref_cids = _sorted_chr_ids(ref_karyo)
    for cid in sorted_ref_cids:
        for gid in ref_karyo[cid]:
            ref_order[_gene_id(gid)] = pos
            pos += 1
    n_ref = pos

    cmap = plt.cm.get_cmap('tab20', max(len(sorted_ref_cids), 1))
    chr_colors = {}
    for i, cid in enumerate(sorted_ref_cids):
        chr_colors[cid] = cmap(i)

    if ref_chrom_of is None:
        ref_chrom_of = {}
        for cid in sorted_ref_cids:
            for gid in ref_karyo[cid]:
                ref_chrom_of[_gene_id(gid)] = cid

    sorted_query_cids = _sorted_chr_ids(query_karyo)
    fig, ax = plt.subplots(figsize=(7, 7))
    q_pos = 0
    xs, ys, colors = [], [], []

    for cid in sorted_query_cids:
        for gid in query_karyo[cid]:
            ref_gene = gene_mapping(_gene_id(gid))
            if ref_gene in ref_order:
                xs.append(ref_order[ref_gene])
                ys.append(q_pos)
                colors.append(chr_colors.get(
                    ref_chrom_of.get(ref_gene, ''),
                    (0.5, 0.5, 0.5, 1.0)))
            q_pos += 1
    n_query = q_pos

    ax.scatter(xs, ys, c=colors, s=5, marker='o',
               alpha=0.6, edgecolors='none')

    ax.set_xlabel(ref_name)
    ax.set_ylabel("{} gene order".format(query_name))
    ax.set_title("{} vs {}".format(query_name, ref_name))
    ax.set_xlim(0, n_ref)
    ax.set_ylim(n_query, 0)

    # Reference chromosome boundaries (vertical)
    ref_labels, _ = _chr_boundaries(ax, sorted_ref_cids, ref_karyo, axis='x')
    # Query chromosome boundaries (horizontal)
    q_labels, _ = _chr_boundaries(ax, sorted_query_cids, query_karyo, axis='y')

    # Chromosome labels on top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([p for p, _ in ref_labels])
    ax2.set_xticklabels([c for _, c in ref_labels],
                        fontsize=7, rotation=45, ha='left')
    ax2.tick_params(length=0)

    # Chromosome labels on right axis
    ax3 = ax.twinx()
    ax3.set_ylim(ax.get_ylim())
    ax3.set_yticks([p for p, _ in q_labels])
    ax3.set_yticklabels([c for _, c in q_labels], fontsize=7)
    ax3.tick_params(length=0)

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def draw_self_dotplot(karyo, outpath, title="Self-dotplot", dpi=200):
    """Draw a self-dotplot for a single karyotype."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        print("WARNING: matplotlib not available, skipping dotplot generation")
        return

    sorted_cids = _sorted_chr_ids(karyo)
    cmap = plt.cm.get_cmap('tab20', max(len(sorted_cids), 1))
    chr_colors = {}
    for i, cid in enumerate(sorted_cids):
        chr_colors[cid] = cmap(i)

    fig, ax = plt.subplots(figsize=(10, 10))
    pos = 0
    xs, ys, colors = [], [], []
    for cid in sorted_cids:
        for gid in karyo[cid]:
            xs.append(pos)
            ys.append(pos)
            colors.append(chr_colors[cid])
            pos += 1
    n_total = pos

    ax.scatter(xs, ys, c=colors, s=6, marker='o',
               alpha=0.6, edgecolors='none')
    ax.set_xlabel("Gene order")
    ax.set_ylabel("Gene order")
    ax.set_title(title)
    ax.set_xlim(0, n_total)
    ax.set_ylim(n_total, 0)

    anc_labels, _ = _chr_boundaries(ax, sorted_cids, karyo, axis='x')
    _chr_boundaries(ax, sorted_cids, karyo, axis='y')

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([p for p, _ in anc_labels])
    ax2.set_xticklabels([c for _, c in anc_labels],
                        fontsize=7, rotation=45, ha='left')
    ax2.tick_params(length=0)

    ax3 = ax.twinx()
    ax3.set_ylim(ax.get_ylim())
    ax3.set_yticks([p for p, _ in anc_labels])
    ax3.set_yticklabels([c for _, c in anc_labels], fontsize=7)
    ax3.tick_params(length=0)

    legend_patches = [Patch(facecolor=chr_colors[c], label=c)
                      for c in sorted_cids]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=7,
              ncol=min(len(sorted_cids), 4))

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def draw_ancestor_dotplots(tree, node_karyotypes, outdir,
                           root_name=None, ancestor_fn=None, dpi=200):
    """Draw dotplots for a root ancestor vs all descendant nodes.

    For each node in the tree (excluding the root), draws a dotplot
    comparing that node's karyotype against the root ancestor's
    karyotype. Also draws a self-dotplot for the ancestor.

    Parameters
    ----------
    tree : ete3.Tree
        Species tree.
    node_karyotypes : dict
        {node_name: {chrom_id: [gene_id, ...]}}
    outdir : str
        Output directory.
    root_name : str, optional
        Name of the root/ancestor node. If None, detected from tree.
    ancestor_fn : callable, optional
        Maps a query gene to its ancestral gene. If None, identity
        is used.
    dpi : int
        Output resolution.
    """
    if root_name is None:
        for node in tree.traverse("preorder"):
            if node.is_root():
                root_name = node.name
                break

    if not root_name or root_name not in node_karyotypes:
        return

    root_karyo = node_karyotypes[root_name]

    # Dotplot for each descendant vs ancestor
    for node in tree.traverse("preorder"):
        name = node.name
        if name == root_name or name not in node_karyotypes:
            continue
        karyo = node_karyotypes[name]
        outpath = os.path.join(outdir, "dotplot_{}.png".format(name))
        draw_dotplot(root_karyo, karyo, outpath,
                     ref_name="Ancestor ({})".format(root_name),
                     query_name=name,
                     gene_mapping=ancestor_fn,
                     dpi=dpi)


def _make_branch_mapping(gene_parent, parent_gene_set):
    """Build a gene mapping function for a parent→child branch dotplot.

    For each child gene, tries to find the corresponding gene in the
    direct parent:
    1. Direct match (same gene ID in parent)
    2. gene_parent chain (WGD/duplicate copies)
    3. Fallback to root ancestor mapping

    Parameters
    ----------
    gene_parent : dict
        {gene_id: immediate_parent_gene_id}
    parent_gene_set : set
        Set of gene IDs present in the parent karyotype.

    Returns
    -------
    callable
        Mapping function: child_gene_id -> parent_gene_id or None.
    """
    def _mapping(gid):
        # 1. Direct match
        if gid in parent_gene_set:
            return gid
        # 2. Follow gene_parent chain
        cur = gid
        visited = set()
        while cur in gene_parent:
            if cur in visited:
                break
            visited.add(cur)
            cur = gene_parent[cur]
            if cur in parent_gene_set:
                return cur
        # 3. No match
        return None
    return _mapping


def draw_branch_dotplots(tree, node_karyotypes, outdir, gene_parent,
                         dpi=200):
    """Draw dotplots for each parent→child branch in the tree.

    For every branch (parent, child), draws a dotplot with the parent
    karyotype as reference (x-axis) and the child karyotype as query
    (y-axis). This reveals the chromosome changes at each evolutionary
    step.

    Parameters
    ----------
    tree : ete3.Tree
        Species tree.
    node_karyotypes : dict
        {node_name: {chrom_id: [gene_id, ...]}}
    outdir : str
        Output directory for branch dotplot PNGs.
    gene_parent : dict
        {gene_id: immediate_parent_gene_id} from GeneTracker.
        Used to map child genes back to parent genes.
    dpi : int
        Output resolution.
    """
    os.makedirs(outdir, exist_ok=True)

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        if parent is None:
            continue
        # Use actual node.name for karyotype lookup, display name for labels
        parent_key = parent.name
        child_key = node.name
        if parent_key not in node_karyotypes or child_key not in node_karyotypes:
            continue
        parent_label = parent_key if parent_key else "root"
        child_label = child_key if child_key else "unnamed"

        parent_karyo = node_karyotypes[parent_key]
        child_karyo = node_karyotypes[child_key]

        # Build set of parent gene IDs for mapping
        parent_gene_set = set()
        for cid in parent_karyo:
            for g in parent_karyo[cid]:
                parent_gene_set.add(_gene_id(g))

        branch_mapping = _make_branch_mapping(gene_parent, parent_gene_set)

        outpath = os.path.join(
            outdir, "dotplot_{}_{}.png".format(parent_label, child_label))
        draw_dotplot(parent_karyo, child_karyo, outpath,
                     ref_name=parent_label,
                     query_name=child_label,
                     gene_mapping=branch_mapping,
                     dpi=dpi)


def draw_sister_dotplots(tree, node_karyotypes, outdir, dpi=200):
    """Draw dotplots between sister nodes (children of the same parent).

    For every internal node with ≥2 children, draws pairwise dotplots
    between all child pairs.

    Parameters
    ----------
    tree : ete3.Tree
        Species tree.
    node_karyotypes : dict
        {node_name: {chrom_id: [gene_id, ...]}}
    outdir : str
        Output directory.
    dpi : int
        Output resolution.
    """
    os.makedirs(outdir, exist_ok=True)

    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        children = [c for c in node.children if c.name in node_karyotypes]
        if len(children) < 2:
            continue
        # Pairwise between sisters
        for i in range(len(children)):
            for j in range(i + 1, len(children)):
                c1_name = children[i].name
                c2_name = children[j].name
                k1 = node_karyotypes[c1_name]
                k2 = node_karyotypes[c2_name]
                outpath = os.path.join(
                    outdir, "dotplot_{}_vs_{}.png".format(c1_name, c2_name))
                draw_dotplot(k1, k2, outpath,
                             ref_name=c1_name,
                             query_name=c2_name,
                             dpi=dpi)


def aag_to_karyo(aag):
    """Convert AncestralAdjacencyGraph to plain karyotype dict.

    Returns
    -------
    dict
        {chrom_id: [gene_id, ...]}
    """
    karyo = {}
    chrom_idx = 0
    for chrom in aag.chromosomes:
        real_nodes = [n for n in chrom if n not in aag.telomeres]
        if not real_nodes:
            continue
        chrom_idx += 1
        cid = "{}_c{}".format(aag.node_id, chrom_idx)
        karyo[cid] = []
        for node in real_nodes:
            if hasattr(node, 'hog_id'):
                karyo[cid].append(node.hog_id)
            else:
                karyo[cid].append(str(node))
    return karyo


def draw_akr_dotplots(akr, outdir, dpi=200):
    """Draw dotplots for AKR reconstructed ancestors vs descendants.

    Draws two sets of plots:
    1. For every internal node, a dotplot against each direct child.
    2. For every node (including leaves), a dotplot against the root.
    3. A self-dotplot for the root ancestor.

    Parameters
    ----------
    akr : AKR
        An AKR instance with ``anc_graphs`` and ``leaf_graphs``
        populated.
    outdir : str
        Output directory.
    dpi : int
        Output resolution.
    """
    os.makedirs(outdir, exist_ok=True)

    # Build karyotype dicts for all nodes
    all_karyos = {}
    for node_id, aag in akr.anc_graphs.items():
        all_karyos[node_id] = aag_to_karyo(aag)
    for node_id, aag in akr.leaf_graphs.items():
        if node_id not in all_karyos:
            all_karyos[node_id] = aag_to_karyo(aag)

    # Find root
    root_name = None
    for node in akr.tree.traverse("preorder"):
        if node.is_root():
            root_name = node.name
            break

    if not root_name or root_name not in all_karyos:
        return

    # Parent-child dotplots
    for node in akr.tree.traverse("preorder"):
        node_id = node.name
        if node_id not in all_karyos:
            continue
        node_karyo = all_karyos[node_id]

        for child in node.children:
            child_id = child.name
            if child_id not in all_karyos:
                continue
            child_karyo = all_karyos[child_id]
            outpath = os.path.join(
                outdir,
                "dotplot_{}_vs_{}.png".format(node_id, child_id))
            draw_dotplot(node_karyo, child_karyo, outpath,
                         ref_name=node_id,
                         query_name=child_id,
                         dpi=dpi)

        # Each non-root node vs root
        if node_id != root_name:
            outpath = os.path.join(
                outdir,
                "dotplot_{}_vs_{}.png".format(root_name, node_id))
            draw_dotplot(all_karyos[root_name], node_karyo, outpath,
                         ref_name=root_name,
                         query_name=node_id,
                         dpi=dpi)
