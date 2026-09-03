# coding: utf-8
"""
paralog: output paralogous gene pairs produced at each branch from HOGs,
         and index collinearity blocks by branch-paralog content.

Internal nodes: sibling child HOGs (same parent) -> genes of the same species
  across different child HOGs are paralogs on that branch.
Leaf species: genes of the same species within the same HOG -> paralogs.
"""
from collections import defaultdict, Counter
import os
from .hog import HOG
from .style import set_species_ticklabels, format_species_label
from . import style
from .mcscan import XCollinearity
from .RunCmdsMP import logger


# ---------------------------------------------------------------------------
#  paralog output (used directly or by indexer)
# ---------------------------------------------------------------------------

class Paralog:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None,
				 nodes=None, species=None, prefix='paralog', **kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.kargs = kargs

	def run(self):
		"""Load HOGs, write paralog TSV, return {branch: frozenset(pairs)}, root_name.

		With hog_tsv: load HOGs from an existing HOGs.tsv instead of rebuilding.
		With inparalog: read an existing --no-index paralog.tsv (first 3 cols).
		"""
		hog_tsv = self.kargs.get('hog_tsv')
		inparalog = self.kargs.get('inparalog')
		fpath = self.prefix + '.paralog.tsv'

		if inparalog:
			branch_pairs, count = self._load_pairs(inparalog)
			root = self._root_name()
			logger.info('Loaded {} paralog pairs from {}'.format(count, inparalog))
			# warn when loaded nodes do not match the species tree nodes
			if self.sptreefile:
				from .tree import number_nodes
				tree_names = {n.name for n in number_nodes(self.sptreefile).traverse()}
				unknown = sorted(set(branch_pairs) - tree_names)
				if unknown:
					logger.warning('Loaded nodes not present in species tree ({}): {} — check that {} is a paralog.tsv (gene1 gene2 node), not a HOGs.tsv'.format(len(unknown), ', '.join(unknown[:10]), inparalog))
			return {b: frozenset(ps) for b, ps in branch_pairs.items()}, root


		self.kargs.update(min_child_species=2, cross_speciation=True, drop_no_cross=True,
						 min_species_ratio=0.2)
		logger.info('reset kargs for paralog: {}'.format(self.kargs))
		hog = HOG(ogfile=self.ogfile, orthfiles=self.orthfiles,
				  sptreefile=self.sptreefile, outpre=self.prefix, **self.kargs)
		if hog_tsv:
			hog.from_tsv(hog_tsv)
		else:
			write_tsv = bool(self.kargs.get('write_hog'))
			if write_tsv:
				hog.outtsv = self.prefix + '.hog.tsv'
			hog.pipe(write_tsv=write_tsv)
		logger.info('Loaded {} HOGs'.format(len(hog.all_hogs)))
		if not hog.all_hogs:
			raise ValueError('0 HOGs after HOG building: check collinearity data and --drop-no-cross (all gene sets dropped)')

		branch_pairs, count = self._write_and_group(hog, fpath)
		logger.info('Output {} paralog pairs to {}'.format(count, fpath))
		root = hog.tree.name if hog.tree else self._root_name()
		return {b: frozenset(ps) for b, ps in branch_pairs.items()}, root

	def _root_name(self):
		"""Root node name from species tree, or 'Root' if no tree."""
		if not self.sptreefile:
			return 'Root'
		from .tree import number_nodes
		return number_nodes(self.sptreefile).name

	def _load_pairs(self, fpath):
		"""Read paralog pairs from an existing paralog.tsv (first 3 columns).

		Applies the same --nodes/--species filters as HOG rebuilding.
		Return ({branch: set of canonical (g1,g2)}, count).
		"""
		node_set = set(self.nodes) if self.nodes else None
		sp_set = set(self.species) if self.species else None
		branch_pairs = defaultdict(set)
		count = 0
		with open(fpath) as f:
			for line in f:
				if line.startswith('#'):
					continue
				parts = line.rstrip().split('	')
				if len(parts) < 3:
					continue
				g1, g2, node = parts[0], parts[1], parts[2]
				if node_set and node not in node_set:
					continue
				if sp_set:
					sp = g1.split('|')[0] if '|' in g1 else g1.rsplit('_', 1)[0]
					if sp not in sp_set:
						continue
				pair = (g1, g2) if g1 < g2 else (g2, g1)
				branch_pairs[node].add(pair)
				count += 1
		return branch_pairs, count

	def _write_and_group(self, hog, fpath):
		"""Write paralog TSV, return {branch: set of canonical (g1,g2)}, count."""
		branch_pairs = defaultdict(set)
		with open(fpath, 'w') as fout:
			fout.write('#gene1\tgene2\tnode\tspecies\tHOG_id\n')
			count = 0
			for g1, g2, node_id, sp, hog_id in hog.iter_branch_paralogs(
					self.nodes, self.species):
				fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
					g1, g2, node_id, sp, hog_id))
				pair = (g1, g2) if g1 < g2 else (g2, g1)
				branch_pairs[node_id].add(pair)
				count += 1
		return branch_pairs, count


# ---------------------------------------------------------------------------
#  branch-paralog indexer (default mode)
# ---------------------------------------------------------------------------

class ParalogIndexer:
	"""Assign collinearity blocks to branches based on paralog content."""

	def __init__(self, ogfile, orthfiles, sptreefile,
				 self_synteny=None, min_n=0, gff=None, min_dist=None,
				 pi_cutoff=0.05, nodes=None, species=None,
				 prefix='paralog_index', box_branch=False, line_plot=False, note=None,
				 tree_heatmap=False, tree_style='distance', heatmap_raster=True,
				 **hog_kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.self_synteny = self_synteny or orthfiles
		self.note = note
		if not self.self_synteny:
			import sys
			logger.error('Paralog indexing needs synteny blocks: '
						 'provide -ss or -s')
			sys.exit(1)
		self.min_n = min_n
		self.gff = gff
		self.min_dist = min_dist
		self.pi_cutoff = pi_cutoff
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.box_branch = box_branch
		self.line_plot = line_plot
		self.tree_heatmap = tree_heatmap
		self.tree_style = tree_style
		self.heatmap_raster = heatmap_raster
		self.hog_kargs = hog_kargs

		# lazy
		self._branch_pairs = None
		self._root_branch = None

	def _load_branch_pairs(self):
		"""Load HOGs and paralog pairs via Paralog."""
		if self._branch_pairs is not None:
			return
		paralog = Paralog(ogfile=self.ogfile, orthfiles=self.orthfiles,
						  sptreefile=self.sptreefile, prefix=self.prefix,
						  nodes=self.nodes, species=self.species,
						  **self.hog_kargs)
		self._branch_pairs, self._root_branch = paralog.run()
		# count paralog pairs per (branch, species)
		self._branch_sp_counts = defaultdict(int)
		for branch, pairs in self._branch_pairs.items():
			for g1, g2 in pairs:
				sp = g1.split('|')[0]
				self._branch_sp_counts[(branch, sp)] += 1

	def _compute_pi(self, block_pairs, branch_pairs):
		"""Return (BPI, n_intersect) for branch paralogy index."""
		if not block_pairs:
			return 0.0, 0
		common = sum(1 for p in block_pairs if p in branch_pairs)
		return common / len(block_pairs), common

	def _canonical_pair(self, g1, g2):
		return (g1, g2) if g1 < g2 else (g2, g1)

	def assign(self):
		"""Assign blocks to branches, return {branch: [(block, sp1, sp2, N, pi), ...]}."""
		self._load_branch_pairs()
		assigned = defaultdict(list)
		threshold = self.pi_cutoff
		root = self._root_branch
		# branch column order: internal nodes first, then leaves, each in
		# tree traversal order; leaves follow tree topology clustering.
		# Without a species tree, fall back to sorted branch order.
		if self.sptreefile:
			from .tree import number_nodes
			tree = number_nodes(self.sptreefile)
			if self.tree_heatmap:
				# Fbranch-style row order: preorder over ALL nodes, so every
				# node's children sit on adjacent rows and the tree panel
				# fans out without crossing connectors
				tree_order = [n.name for n in tree.traverse('preorder')]
			else:
				internal = [n.name for n in tree.traverse() if not n.is_leaf()]
				leaves = tree.get_leaf_names()[::-1]  # reverse: outgroups first
				tree_order = internal + leaves
			# 列所有物种树 branch（含没旁系的，全 0 也占一列，别缺列）；树外 branch 兜底
			branches = list(tree_order)
			branches += sorted(set(self._branch_pairs) - set(branches))
			# --species 时只保留该物种相关的 branch（该物种 leaf + 含该物种的内部节点）
			if self.species:
				sp_set = set(self.species)
				keep = []
				for b in branches:
					node = next((n for n in tree.traverse() if n.name == b), None)
					if node is None:
						continue
					if node.is_leaf():
						if node.name in sp_set:
							keep.append(b)
					else:
						# internal: keep if any descendant leaf species in sp_set
						if sp_set & set(node.get_leaf_names()):
							keep.append(b)
				branches = keep
		else:
			branches = sorted(self._branch_pairs)
		self._pi_branches = branches
		self._ordered_branches = list(branches)  # internal+leaves in tree order; shared by heatmap/stats
		self._pi_rows = []  # (block_id, N, pi_vector) for heatmap
		self._branch_raw_signal = defaultdict(int)  # raw paralog count per branch (all blocks, even weak)

		# tree leaf species for consistency check (assign must never produce
		# stats keys for species absent from the tree — write_stats sorts by
		# _ordered_branches.index(species) and would raise ValueError)
		tree_species = None
		if self.sptreefile:
			from .tree import number_nodes
			tree_species = set(number_nodes(self.sptreefile).get_leaf_names())
		excluded_species = defaultdict(lambda: [0, 0])  # sp -> [n_blocks, n_genes]

		for rc in XCollinearity(self.self_synteny, gff=self.gff):
			if rc.N < self.min_n:
				continue
			if rc.species1 != rc.species2:
				continue
			if self.min_dist and rc.is_tandem(self.min_dist):
				continue
			if self.species and rc.species1 not in self.species:
				continue
			if tree_species is not None and rc.species1 not in tree_species:
				# single warning after the loop; never touch filenames
				ex = excluded_species[rc.species1]
				ex[0] += 1
				ex[1] += rc.N
				continue

			block_pairs = [self._canonical_pair(g1, g2) for g1, g2 in rc.pairs]

			best_branch = root
			best_pi = 0.0
			best_nparalog = 0
			pi_vector = []
			for branch in branches:
				pi, n_paralog = self._compute_pi(block_pairs,
												 self._branch_pairs.get(branch, set()))
				pi_vector.append(pi)
				self._branch_raw_signal[branch] += n_paralog
				if pi > best_pi:
					best_pi = pi
					best_nparalog = n_paralog
					best_branch = branch

			if best_pi < threshold:
				best_branch = "Unassigned"

			# store data immediately — rc is a shared mutable object
			assigned[best_branch].append(
				(rc.block, rc.species1, rc.N, best_nparalog, best_pi))
			self._pi_rows.append((rc.id, rc.species1, rc.N,
								  best_branch, tuple(pi_vector)))
		self._pi_branches = branches
		if excluded_species:
			logger.warning('Excluded {} self-synteny blocks ({} genes) from {} species not in species tree: {}'.format(
				sum(v[0] for v in excluded_species.values()),
				sum(v[1] for v in excluded_species.values()),
				len(excluded_species),
				', '.join('{} ({} blocks)'.format(sp, v[0]) for sp, v in sorted(excluded_species.items()))))
		logger.info('Assigned blocks to {} branches'.format(len(assigned)))
		n_assigned_blocks = sum(len(v) for b, v in assigned.items() if b != 'Unassigned')
		n_unassigned_blocks = len(assigned.get('Unassigned', []))
		n_assigned_genes = sum(it[2] for b, v in assigned.items() if b != 'Unassigned' for it in v)
		n_unassigned_genes = sum(it[2] for it in assigned.get('Unassigned', []))
		logger.info('Blocks assigned/unassigned: {}/{}; genes assigned/unassigned: {}/{}'.format(
			n_assigned_blocks, n_unassigned_blocks, n_assigned_genes, n_unassigned_genes))
		if assigned and n_assigned_blocks == 0:
			logger.warning('All {} blocks unassigned (no branch reached PI >= {}): check that the paralog input matches the species tree and collinearity data'.format(
				n_unassigned_blocks, threshold))
		return assigned

	def write_heatmap_matrix(self, assigned):
		"""Write block x branch PI matrix TSV (always, even without heatmap plot).

		Rows (blocks) are sorted by their PI vector (descending) so blocks
		with similar branch-PI profiles end up adjacent.  With cluster=True
		rows are ordered by hierarchical clustering instead.
		Column 'assigned_branch' records which branch each block was
		assigned to.
		"""
		fpath = self.prefix + '.heatmap.tsv'
		branches = self._ordered_branches
		branch_idx = {b: i for i, b in enumerate(branches)}
		if getattr(self, 'heatmap_cluster', False):
			# hierarchical clustering of rows
			import numpy as np
			from scipy.cluster.hierarchy import linkage, leaves_list
			M = np.array([r[4] for r in self._pi_rows], dtype=float)
			if M.size:
				Z = linkage(M, method='ward')
				order = leaves_list(Z)
				rows = [self._pi_rows[i] for i in order]
			else:
				rows = []
		else:
			# sort by assigned-branch column index, then PI vector (desc);
			# Unassigned blocks (no paralog signal) go first
			rows = sorted(self._pi_rows,
						  key=lambda r: (0 if r[3] == "Unassigned" else 1,
										 branch_idx.get(r[3], len(branches)),
										 tuple(-v for v in r[4])))

		# TSV output: ordered by (branch, species, block), Unassigned last;
		# independent of the figure row order (rows/_heatmap_rows keep figure order)
		tsv_rows = sorted(rows,
					  key=lambda r: (branch_idx.get(r[3], len(branches)),
									 self._ordered_branches.index(r[1]),
									 r[0]))
		with open(fpath, 'w') as fout:
			fout.write('#block\tspecies\tN\tassigned_branch\t' + '\t'.join(branches) + '\n')
			for bid, sp, N, bch, vec in tsv_rows:
				fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
					bid, sp, N, bch, '\t'.join('{:.4f}'.format(v) for v in vec)))
		logger.info('Heatmap matrix written to {} ({} blocks x {} branches)'.format(
			fpath, len(rows), len(branches)))
		self._heatmap_rows = rows
		return rows

	def write_heatmap(self, assigned):
		"""Draw the heatmap figure from the PI matrix."""
		fpath = self.prefix + '.heatmap.tsv'
		branches = self._ordered_branches
		rows = self._heatmap_rows
		branch_idx = {b: i for i, b in enumerate(branches)}


		# figure
		import numpy as np
		import matplotlib.pyplot as plt
		M = np.array([r[4] for r in rows], dtype=float)
		if M.size == 0:
			logger.warning('No blocks for heatmap')
			return

		# column widths: uniform, or scaled by gene count when --scale given
		# (transposed layout: columns are blocks)
		Ns = np.array([r[2] for r in rows], dtype=float)
		method = getattr(self, 'heatmap_scale', None)
		if method:
			if method == 'linear':
				w = Ns
			elif method == 'log2':
				w = np.log2(1 + Ns)
			elif method == 'log10':
				w = np.log10(1 + Ns)
			elif method == 'sqrt':
				w = np.sqrt(Ns)
			else:  # 'log' = natural log
				w = np.log1p(Ns)
			w = w / w.max() #* (1 - 1e-3) + 1e-3  # tiny floor
		else:
			w = np.ones(len(rows))
		# cumulative x positions (left of each column)
		x = np.concatenate([[0.0], np.cumsum(w)])

		fig_w = min(7, max(4, 0.02 * M.shape[0])) + (1.8 if self.tree_heatmap else 0)
		fig_h = min(7, max(3, 0.4 * len(branches)))
		fig = plt.figure(figsize=(fig_w, fig_h))
		cmap = plt.get_cmap('YlOrRd')
		from matplotlib import gridspec
		# second row (cbar + legend) needs a fixed physical height regardless
		# of branch count, else short figures squeeze it onto the first row
		cbar_need = 1.1  # inches: legend two rows + cbar + xlabels clearance
		hr = max(1.0, len(branches) * cbar_need / max(fig_h - cbar_need, 0.1))
		# tree panel on the left column; heatmap columns shift right by one
		ncol_tree = 1 if self.tree_heatmap else 0
		if ncol_tree:
			if self.line_plot:
				gs = gridspec.GridSpec(2, 3, width_ratios=[2, 6, 2],
										height_ratios=[len(branches), hr],
										wspace=0.05, hspace=0.05)
				ax_tree = fig.add_subplot(gs[0, 0])
				ax = fig.add_subplot(gs[0, 1])
				ax_line = fig.add_subplot(gs[0, 2])
				ax_leg = fig.add_subplot(gs[1, 2])
				cax = fig.add_subplot(gs[1, 1])
			else:
				gs = gridspec.GridSpec(2, 2, width_ratios=[2, 6],
										height_ratios=[len(branches), hr],
										wspace=0.05, hspace=0.05)
				ax_tree = fig.add_subplot(gs[0, 0])
				ax = fig.add_subplot(gs[0, 1])
				ax_line = None
				ax_leg = None
				cax = fig.add_subplot(gs[1, 1])
		else:
			if self.line_plot:
				gs = gridspec.GridSpec(2, 2, width_ratios=[6, 2],
										height_ratios=[len(branches), hr],
										wspace=0.05, hspace=0.05)
				ax = fig.add_subplot(gs[0, 0])
				ax_line = fig.add_subplot(gs[0, 1])
				ax_leg = fig.add_subplot(gs[1, 1])
				cax = fig.add_subplot(gs[1, 0])
			else:
				gs = gridspec.GridSpec(2, 1, height_ratios=[len(branches), hr])
				ax = fig.add_subplot(gs[0, 0])
				cax = fig.add_subplot(gs[1, 0])
				ax_line = None
		# transposed: one row per branch, segments = blocks
		for i, branch_vec in enumerate(M.T):
			colors = [cmap(v) for v in branch_vec]
			ax.broken_barh(list(zip(x[:-1], w)), (i - 0.5, 1.0),
						   facecolors=colors, edgecolors='none',
						   linewidths=0, rasterized=self.heatmap_raster)
		ax.set_xlim(0, x[-1])
		ax.set_ylim(-0.5, len(branches) - 0.5)
		ax.invert_yaxis()  # branches top-down (first branch at top)
		ax.set_yticks(range(len(branches)))
		tick_fs = max(6, min(12, 240 // max(len(branches), 1)))
		# no y label when the tree panel sits on the left
		if not self.tree_heatmap:
			ax.set_ylabel('Branch', fontsize=14)
		if ax_line is not None:
			ax.set_yticklabels([])
			ax.tick_params(axis='y', left=False, labelleft=False)
			ax_line.yaxis.tick_right()  # branch names on the far right
			set_species_ticklabels(ax_line, branches, fontsize=tick_fs)
		else:
			ax.yaxis.tick_right()  # branch names on the right
			set_species_ticklabels(ax, branches, fontsize=tick_fs)
			ax.tick_params(axis='y', pad=1.75)  # half default(3.5)
		ax.set_xticks([])
		ax.set_xlabel('Syntenic block', fontsize=14)
		# optional: grey dashed rectangles around blocks assigned to each branch
		if self.box_branch:
			from collections import defaultdict
			from matplotlib.patches import Rectangle
			blocks_by_branch = defaultdict(list)
			for j, (bid, sp, N, bch, vec) in enumerate(rows):
				blocks_by_branch[bch].append(j)
			for bch, cols in blocks_by_branch.items():
				i = branch_idx.get(bch)
				if i is None:
					continue
				x0 = x[cols[0]]
				x1 = x[cols[-1] + 1]
				rect = Rectangle((x0, i - 0.5), x1 - x0, 1.0,
								 fill=False, edgecolor='0.5', lw=0.6,
								 linestyle=':', zorder=10)
				ax.add_patch(rect)
		# optional: line plot on the right showing raw and assigned paralog signal per branch
		if ax_line is not None:
			raw_sig = [self._branch_raw_signal.get(b, 0) for b in branches]
			asgn_sig = [sum(it[2] for it in assigned.get(b, [])) if b in assigned else 0
						for b in branches]
			ax_line.plot(raw_sig, range(len(branches)), 'b-', lw=1.2, alpha=0.8,
						 label='Raw paralog pairs')
			ax_line.plot(asgn_sig, range(len(branches)), 'orange', lw=1.2, alpha=0.8,
						 label='Assigned gene pairs')
			ax_line.set_ylim(-0.5, len(branches) - 0.5)
			ax_line.invert_yaxis()  # match heatmap direction
			ax_line.set_xlabel('Number of gene pairs', fontsize=9)
			ax_line.tick_params(axis='y', labelleft=False, left=False, pad=6)  # half of 12: label closer to line panel
			# ax_line.set_xticks([min(raw_sig + asgn_sig), max(raw_sig + asgn_sig)])
			from matplotlib.ticker import AutoMinorLocator
			ax_line.xaxis.set_minor_locator(AutoMinorLocator())
			ax_leg.axis('off')
			lg = ax_leg.legend(ax_line.get_lines(), ['Inferred paralog pairs', 'Assigned syntenic gene pairs'],
							 fontsize=8, loc='center', ncol=1)
		import matplotlib as mpl
		fig.tight_layout()
		# tree panel: draw after tight_layout so it uses final axes geometry
		if ncol_tree:
			if self.sptreefile and os.path.exists(self.sptreefile):
				self._draw_tree_panel(ax_tree, branches)
			else:
				logger.warning('--tree-heatmap requires a valid -t tree file; panel left empty')
		# cbar top must clear the first row's tight bbox (xlabel included),
		# else it lands on the 'Syntenic block'/'Number of gene pairs' labels
		fig.canvas.draw()
		renderer = fig.canvas.get_renderer()
		floor = min(ax.get_tightbbox(renderer).y0 for ax in (ax, ax_line) if ax is not None)
		floor = floor / fig.bbox.height  # display -> figure coords
		fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1),
										   cmap=cmap),
					 cax=cax, label='BPI', orientation='horizontal')
		cbar_h = 0.02
		cbar_top = cax.get_position().y1 - 0.02  # current anchor from gridspec
		if ax_line is not None:
			lg_bbox = lg.get_window_extent().transformed(fig.transFigure.inverted())
			# the legend itself must also clear the first row; on short
			# figures its centered position in ax_leg can poke above the floor
			pos = ax_leg.get_position()
			if lg_bbox.y1 > floor - 0.02:
				dy = lg_bbox.y1 - (floor - 0.02)
				ax_leg.set_position([pos.x0, pos.y0 - dy, pos.width, pos.height])
				fig.canvas.draw()
				lg_bbox = lg.get_window_extent().transformed(fig.transFigure.inverted())
			cbar_top = min(cbar_top, lg_bbox.y1)
		cbar_top = min(cbar_top, floor - 0.02)  # never overlap the first row
		cax.set_position([cax.get_position().x0, cbar_top - cbar_h,
						  0.4 * cax.get_position().width, cbar_h])
		if getattr(self, 'note', None):
			fig.text(0.1, 0.1, self.note, fontsize=8, va='bottom', ha='left',
					 color='0.4', bbox=dict(facecolor='white', alpha=1, edgecolor='none', pad=2))
		fig.savefig(self.prefix + '.heatmap.pdf', dpi=300, bbox_inches='tight')
		fig.savefig(self.prefix + '.heatmap.png', dpi=150, bbox_inches='tight')
		plt.close(fig)
		logger.info('Heatmap written to {}.heatmap.pdf/png'.format(self.prefix))

	def write_stats(self, assigned):
		"""Write per-branch per-species statistics TSV, return stats dict."""
		fpath = self.prefix + '.paralog.stats.tsv'
		stats = defaultdict(lambda: [0, 0, 0, 0.0])
		for branch, items in assigned.items():
			for block, sp, N, n_paralog, pi in items:
				if sp is None:
					continue
				s = stats[(branch, sp)]
				s[0] += 1
				s[1] += N
				s[2] += n_paralog
				s[3] += pi

		with open(fpath, 'w') as fout:
			fout.write('#branch	species	assigned_paralog_pairs	assigned_syntenic_blocks	'
								   'gene_pairs_in_syntenic_blocks	paralog_pairs_in_syntenic_blocks	'
								   'mean_PI	weighted_PI\n')
			ordered = {b: i for i, b in enumerate(self._ordered_branches)}
			for (branch, sp), (blocks, gp, pp, sum_pi) in sorted(
					stats.items(),
					key=lambda kv: (ordered.get(kv[0][0], len(ordered)),  # Unassigned last
								  self._ordered_branches.index(kv[0][1]))):
				mean_pi = sum_pi / blocks if blocks else 0.0
				n_paralogs = self._branch_sp_counts.get((branch, sp), 0)
				wpi = pp / gp if gp else 0.0
				fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\n'.format(
					branch, sp, n_paralogs, blocks, gp, pp, mean_pi, wpi))
		logger.info('Stats written to {}'.format(fpath))
		return stats

	def _draw_tree_panel(self, ax, branches):
		"""Draw the species tree with rows aligned to heatmap branch rows.

		Fbranch-style: every branch (internal + leaf) sits on its own heatmap
		row; dotted guides extend internal-node rows to the panel edge.
		"""
		from .tree import number_nodes, TreePlotter
		tree = number_nodes(self.sptreefile)
		tp = TreePlotter(tree, row_names=branches,
						 internal_node_to_present=True,
						 tree_style=self.tree_style,
						 leaf_dotted_color='tab:green')
		tp.plot(ax)

	def plot_tree(self, assigned, stats):
		"""Render species tree with pies: [syntenic, non-syntenic] paralog pairs.

		Pie size maps to total paralog_pairs per branch; pie fractions show
		the proportion of paralog pairs found in syntenic blocks.
		"""
		import os
		if 'QT_QPA_PLATFORM' not in os.environ:
			os.environ['QT_QPA_PLATFORM'] = 'offscreen'
		from ete3 import TreeStyle, PieChartFace, TextFace, AttrFace, NodeStyle, faces as etefaces
		from .tree import number_nodes

		# aggregate per branch over species
		branch_syn = defaultdict(int)   # branch -> syntenic paralog pairs
		branch_all = defaultdict(int)   # branch -> total paralog pairs
		for (branch, sp), (blocks, gp, pp, sum_pi) in stats.items():
			branch_syn[branch] += pp
			branch_all[branch] += self._branch_sp_counts.get((branch, sp), 0)

		tree = number_nodes(self.sptreefile)
		# number of descendant species per node — normalize pie size
		node_nspecies = {}
		for node in tree.traverse():
			node_nspecies[node.name] = len(node.get_leaf_names())

		ns = NodeStyle()
		ns['size'] = 0
		ns['hz_line_width'] = 1
		ns['vt_line_width'] = 1

		# scale pie size by paralog pairs per descendant species
		branch_per_sp = {b: 1.0 * n / max(node_nspecies.get(b, 1), 1)
						 for b, n in branch_all.items()}
		max_n = max(branch_per_sp.values()) if branch_per_sp else 1
		scale = 60.0 / max_n  # max pie 60px

		def layout(node):
			node.set_style(ns)
			nid = node.name
			if nid in branch_all:
				all_n = branch_all[nid]
				if all_n == 0:
					return
				syn_n = branch_syn.get(nid, 0)
				non_n = all_n - syn_n
				pcts = [100.0 * syn_n / all_n, 100.0 * non_n / all_n]
				size = max(8, int(scale * branch_per_sp[nid]))  # pie diameter
				pie = PieChartFace(pcts, size, size,
								   colors=['#377eb8', '#d9d9d9'],
								   line_color=None)
				pie.opacity = 0.9
				if node.is_leaf():
					etefaces.add_face_to_node(pie, node, column=0)
				else:
					etefaces.add_face_to_node(pie, node, column=0, position='branch-right')
			if node.is_leaf():
				text, st = format_species_label(node.name)
				name_face = TextFace(text, fsize=10, fgcolor='black',
									 ftype=style.DEFAULT_FONT, fstyle=st)
				etefaces.add_face_to_node(name_face, node, column=1)
			else:
				name_face = TextFace(nid, fsize=7, fgcolor='#888888', ftype=style.DEFAULT_FONT)
				etefaces.add_face_to_node(name_face, node, column=0, position='branch-top')

		ts = TreeStyle()
		ts.layout_fn = layout
		ts.show_leaf_name = False
		ts.scale = 200
		ts.branch_vertical_margin = 10

		out_pdf = self.prefix + '.tree.pdf'
		out_png = self.prefix + '.tree.png'
		tree.render(out_png, tree_style=ts)
		tree.render(out_pdf, tree_style=ts)
		logger.info('Tree plot written to {}'.format(out_pdf))

	def write_blocks(self, assigned):
		"""Write assigned blocks per branch."""
		for branch, items in assigned.items():
			fpath = '{}.{}.blocks'.format(self.prefix, branch)
			with open(fpath, 'w') as fout:
				for block, sp, N, n_paralog, pi in items:
					fout.write(block)
		logger.info('Block files written for {} branches'.format(len(assigned)))




# ---------------------------------------------------------------------------
#  CLI entry
# ---------------------------------------------------------------------------

def xmain(**kargs):
	no_index = kargs.pop('no_index', False)
	hog_tsv = kargs.get('hog_tsv')
	inparalog = kargs.get('inparalog')
	if not no_index and not hog_tsv and not inparalog:
		if not kargs.get('ogfile') or not kargs.get('orthfiles'):
			import sys
			logger.error('Either -og/-s or --hog/--paralog must be provided')
			sys.exit(1)
	if no_index:
		kargs.pop('output', None)
		Paralog(**kargs).run()
	else:
		tree_plot = kargs.pop('tree_plot', False)
		heatmap = kargs.pop('heatmap', False)
		heatmap_cluster = kargs.pop('heatmap_cluster', False)
		heatmap_scale = kargs.pop('heatmap_scale', None)
		heatmap_raster = kargs.pop('heatmap_raster', True)
		indexer = ParalogIndexer(**kargs)
		indexer.heatmap_cluster = heatmap_cluster
		indexer.heatmap_scale = heatmap_scale
		indexer.heatmap_raster = heatmap_raster
		assigned = indexer.assign()
		stats = indexer.write_stats(assigned)
		indexer.write_heatmap_matrix(assigned)
		indexer.write_blocks(assigned)
		if tree_plot:
			if kargs.get('sptreefile'):
				indexer.plot_tree(assigned, stats)
			else:
				logger.warning('--tree-plot requires -t; skipped')
		if heatmap:
			indexer.write_heatmap(assigned)
