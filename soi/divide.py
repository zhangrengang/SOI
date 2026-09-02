# -*- coding: utf-8 -*-
"""DivideGraph refine method (GETHOGs 1, Altenhoff et al. 2013).

Faithful re-implementation of Algorithm 3/4 from:
	Altenhoff AM, Gil M, Gonnet GH, Dessimoz C. "Inferring hierarchical
	orthologous groups from orthologous gene pairs." PLoS ONE 8(1): e53786 (2013).

Recursive minimum-cut splitting of an orthology graph:
	- FRITS (Fraction Reachable In Two Steps, Algorithm 4) decides WHETHER to
	  cut: cut only when FRITS < alpha (default 0.6).  FRITS is weight-blind:
	  it only looks at graph topology.
	- When triggered, a weighted minimum cut (Stoer-Wagner, deterministic)
	  decides WHERE to cut.  Input weights are used as-is (OI values are
	  already normalized 0-1; the spec's sp_ratio/sp_z normalization targets
	  bit-score scales and is intentionally not implemented here).
	- Diagnostic flags on every cut (spec section 5, our extension):
		near_tie   : cut_cost < gamma * w_iso -> split is barely cheaper than
		             isolating the cheapest single node; evidence weak.
		frits_pass : FRITS >= alpha on a graph that still looks bimodal.

Implementation notes (deviations from the paper, per spec section 3):
	- Stoer-Wagner instead of Karger-Stein (deterministic, reproducible;
	  networkx requires a connected graph, hence the connected-component
	  decomposition in 3.5 which the paper also performs).
	- Degenerate cuts (one side empty) are treated as un-splittable (3.4).
	- This module exposes divide_components(G, ...) with the same return
	  shape as soi.hog.mcl_components: list of node sets.  It plugs into the
	  refine_fn slot of the HOG splitting loop.

Known failure modes on WGD-derived graphs (spec 4.5/4.6, verified on the
monocots test case 2026-09-02):
	- sparse-bridge merges (diameter ~2) never trigger the cut: FRITS high
	  -> under-splitting (the "small-graph blindness" boundary);
	- when the cut IS triggered on dense blocks, the parsimony minimum
	  prefers isolating a single node (cost O(1)) over separating copy
	  lineages (cost O(species)); singleton over-splitting.
	The near_tie flag surfaces exactly these ambiguous cuts.
"""

import networkx as nx

from .RunCmdsMP import logger


def frits(G, sample=100, seed=None):
	"""Fraction of nodes reachable within 2 steps, averaged over sampled nodes.

	Faithful to Algorithm 4.  Weight-blind.  sample=None means all nodes
	(exhaustive); otherwise a random subsample (spec uses constant c=100).
	"""
	import random
	nodes = list(G.nodes())
	n = len(nodes)
	if n == 0:
		return 1.0
	if sample is not None and sample < n:
		rng = random.Random(seed)
		starts = rng.sample(nodes, sample)
	else:
		starts = nodes
	adj = {v: set(G[v]) for v in nodes}
	total = 0
	for v in starts:
		reached = {v}
		frontier = {v}
		for _ in range(2):
			new = set()
			for w in frontier:
				new |= adj[w] - reached
			if not new:
				break
			frontier = new
			reached |= new
		total += len(reached)
	return total / (len(starts) * n)


def w_iso(G):
	"""Cheapest 'isolate one node' alternative: min over nodes of the sum of
	that node's incident edge weights (spec section 5)."""
	cheapest = None
	for v in G.nodes():
		s = sum(d.get('weight', 1.0) for _, _, d in G.edges(v, data=True))
		if cheapest is None or s < cheapest:
			cheapest = s
	return cheapest if cheapest is not None else 0.0


def _mincut(G):
	"""Weighted minimum cut on a CONNECTED graph.  Returns (cost, setA, setB).

	Raises if the graph is not connected (caller decomposes first)."""
	cost, (a, b) = nx.stoer_wagner(G, weight='weight')
	return cost, set(a), set(b)


def divide_components(G, alpha=0.6, gamma=1.2, sample=100, seed=None,
					  diagnostics=None, _depth=0, _component=True):
	"""Refine a graph into subgraphs via recursive minimum cuts.

	Same contract as mcl_components: input nx.Graph with 'weight' edge
	attributes, output list of node sets (singletons included).

	diagnostics: optional list; every cut appends a dict record
	(component_id, n_nodes, n_edges, reachability_before, alpha, cut_cost,
	cut_edges, w_iso, flag, depth).
	"""
	nodes = list(G.nodes())
	if not nodes:
		return []
	if len(nodes) == 1:
		return [set(nodes)]

	# spec 3.5: decompose into connected components first (Stoer-Wagner
	# needs connectivity; the paper also never merges across components)
	if _component:
		out = []
		for cc in nx.connected_components(G):
			out.extend(divide_components(G.subgraph(cc).copy(), alpha=alpha,
										 gamma=gamma, sample=sample, seed=seed,
										 diagnostics=diagnostics,
										 _depth=_depth, _component=False))
		return out

	r = frits(G, sample=sample, seed=seed)
	if r >= alpha:
		# spec: conforming graph -> stop.  frits_pass hints possible
		# under-splitting on a bimodal graph (our diagnostic, not the paper's)
		if diagnostics is not None and _depth == 0:
			diagnostics.append({'n_nodes': len(nodes),
								'n_edges': G.number_of_edges(),
								'reachability_before': r, 'alpha': alpha,
								'flag': 'frits_pass', 'cut_cost': None,
								'cut_edges': None, 'w_iso': None,
								'depth': _depth})
		return [set(nodes)]

	try:
		cost, a, b = _mincut(G)
	except (nx.NetworkXError, ValueError) as e:
		logger.debug('DivideGraph: mincut failed on %d nodes: %s' % (len(nodes), e))
		return [set(nodes)]
	if not a or not b:
		# spec 3.4: degenerate cut -> un-splittable
		return [set(nodes)]

	diag = {'n_nodes': len(nodes), 'n_edges': G.number_of_edges(),
			'reachability_before': r, 'alpha': alpha, 'cut_cost': cost,
			'cut_edges': sorted((u, v, d.get('weight', 1.0))
								for u, v, d in G.edges(data=True)
								if (u in a) != (v in a)),
			'w_iso': w_iso(G), 'depth': _depth}
	diag['flag'] = 'near_tie' if cost < gamma * diag['w_iso'] else None
	if diagnostics is not None:
		diagnostics.append(diag)
	if diag['flag']:
		logger.debug('DivideGraph cut flagged %s: cost=%.3f w_iso=%.3f (%d nodes)'
					 % (diag['flag'], cost, diag['w_iso'], len(nodes)))

	out = []
	for sub in (a, b):
		out.extend(divide_components(G.subgraph(sub).copy(), alpha=alpha,
									 gamma=gamma, sample=sample, seed=seed,
									 diagnostics=diagnostics,
									 _depth=_depth + 1, _component=False))
	return out
