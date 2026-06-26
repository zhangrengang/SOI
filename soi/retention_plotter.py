# coding: utf-8
"""retention: gene retention rate along chromosomes.

Computes the proportion of reference genes that have syntenic orthologs
in each query species, using sliding-window and per-chromosome line plots
with multiple query species overlaid in different colours.
"""
import sys
import argparse
from collections import OrderedDict, Counter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from .mcscan import XCollinearity, XGff
from .small_tools import sort_version
from .RunCmdsMP import logger

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42

# Colour palette for multiple queries
_PALETTE = ['#2166ac', '#b2182b', '#4daf4a', '#ff7f00', '#984ea3',
			'#a65628', '#f781bf', '#999999', '#e41a1c', '#377eb8']


def retention_args(parser):
	parser.add_argument('-s', '-synteny', metavar='FILE', type=str, required=True, nargs='+',
						dest='collinearity',
						help="Syntenic block file (*.collinearity, output of MCSCANX/WGDI)[required]")
	parser.add_argument('-g', '-gff', metavar='FILE', type=str, required=True, nargs='+',
						dest='gff',
						help="Gene annotation gff file (*.gff, one of MCSCANX/WGDI input)[required]")
	parser.add_argument('-r', '-ref', metavar='SPECIES', type=str, required=True, dest='ref',
						help="Reference species [required]")
	parser.add_argument('-q', '-qry', metavar='SPECIES', nargs='+', type=str, dest='qry',
						required=True, help="Query species [required]")
	parser.add_argument('-pre', '-prefix', metavar='PREFIX', type=str,
						dest='output', default=None, help="Output prefix")
	parser.add_argument('--format', metavar='FMT', action='append',
						default=['pdf', 'png'], help="Figure output format [default=%(default)s]")
	parser.add_argument('--window_size', metavar='INT', type=int, default=50,
						help="Window size (number of genes) [default=%(default)s]")
	parser.add_argument('--window_step', metavar='INT', type=int, default=1,
						help="Window step (number of genes) [default=%(default)s]")
	parser.add_argument('--include', metavar='FILE', type=str, default=None,
						help="Gene list for denominator (only count these genes)")
	parser.add_argument('--exclude', metavar='FILE', type=str, default=None,
						help="Gene list to exclude from denominator")
	parser.add_argument('--chrs', metavar='CHR', nargs='+', type=str, default=None,
						help="Chromosomes to plot (default: all with >= min_genes)")
	parser.add_argument('--min_genes', metavar='INT', type=int, default=200,
						help="Minimum genes per chromosome to include [default=%(default)s]")
	parser.add_argument('--count-duplicates', action='store_true', default=False,
						dest='count_duplicates',
						help="Count syntenic genes with multiplicity "
							 "(for WGD query, retention may exceed 1)")
	parser.add_argument('--min_same_block', type=int, default=25,
						help=argparse.SUPPRESS)


class Args:
	def __init__(self):
		pass


def xmain(**kargs):
	args = Args()
	for k, v in kargs.items():
		setattr(args, k, v)
	return main(args)


def main(args):
	if args.output is None:
		sps = [args.ref] + args.qry
		short = [x[:2] for x in sps]
		args.output = 'ret.' + '-'.join(short) + '_w' + str(args.window_size)
	format = args.format if isinstance(args.format, list) else [args.format]
	logger.info('{} x {} species'.format(len(args.qry), args.ref))

	# Load synteny -> syntenic gene counts per query
	syn_genes = parse_collinearity(args.collinearity, args.ref, args.qry,
									min_same_block=getattr(args, 'min_same_block', 25))
	if not args.count_duplicates:
		syn_sets = {}
		for qry, counter in syn_genes.items():
			syn_sets[qry] = set(counter.keys())
		syn_data = syn_sets
	else:
		syn_data = syn_genes

	# Load GFF
	d_chrom_genes, d_chrom_paths = parse_gff(args.gff, args.ref)

	# Optional gene filtering
	include_set = _load_gene_list(args.include) if args.include else None
	exclude_set = _load_gene_list(args.exclude) if args.exclude else None

	# Determine chromosomes
	if args.chrs:
		# User-specified order, no sorting
		plot_chroms = [c for c in args.chrs if c in d_chrom_paths]
		if not plot_chroms:
			logger.error('None of the specified chromosomes found in GFF')
			sys.exit(1)
	else:
		candidates = [c for c in d_chrom_paths if len(d_chrom_paths[c]) >= args.min_genes]
		plot_chroms = sort_version(candidates)

	logger.info('Chromosomes to plot: {} ({} total, min_genes={})'.format(
		len(plot_chroms), len(d_chrom_paths), args.min_genes))

	# Sliding-window retention: qry -> {chrom -> [(pos, rate), ...]}
	all_data = OrderedDict()
	for qry in args.qry:
		syn = syn_data.get(qry, set() if not args.count_duplicates else Counter())
		chrom_data = OrderedDict()
		for chrom in plot_chroms:
			path = d_chrom_paths[chrom]
			points = sliding_retention(path, syn, args.window_size,
										args.window_step, include_set, exclude_set,
										count_duplicates=args.count_duplicates)
			if points:
				chrom_data[chrom] = points
		all_data[qry] = chrom_data

	# Plot: one row per chromosome, one column; all query lines in same panel
	n_chroms = len(plot_chroms)
	figsize = (6, 1.2 * n_chroms)
	fig, axes = plt.subplots(n_chroms, 1, figsize=figsize,
							  sharex=True, sharey=True, squeeze=False)

	# Global y-max and per-chromosome x-max for legend placement
	ymax = 0
	chrom_xmax = {}  # chrom -> max gene index
	for qry in args.qry:
		for chrom in plot_chroms:
			pts = all_data.get(qry, {}).get(chrom, [])
			if pts:
				ymax = max(ymax, max(r for _, r in pts))
				chrom_xmax[chrom] = max(chrom_xmax.get(chrom, 0), max(p for p, _ in pts))
	if not args.count_duplicates:
		ymax = min(1.0, ymax * 1.05) if ymax > 0 else 1.0
	else:
		ymax = ymax * 1.05 if ymax > 0 else 1.0

	# Put legend on the chromosome with the shortest x-range (most empty space)
	if chrom_xmax:
		legend_chrom = min(chrom_xmax, key=lambda c: chrom_xmax[c])
	else:
		legend_chrom = plot_chroms[0] if plot_chroms else None

	for ci, chrom in enumerate(plot_chroms):
		ax = axes[ci][0]
		for qi, qry in enumerate(args.qry):
			pts = all_data.get(qry, {}).get(chrom, [])
			if not pts:
				continue
			positions = [p for p, _ in pts]
			rates = [r for _, r in pts]
			color = _PALETTE[qi % len(_PALETTE)]
			ax.step(positions, rates, color=color, linewidth=0.8, label=qry, where='mid')
			ax.axhline(y=np.mean(rates), color=color, linewidth=0.5,
						linestyle='--', alpha=0.6)
		ax.set_ylim(0, ymax)
		ax.set_title(chrom, fontsize=9)
		ax.set_ylabel('Retention', fontsize=8)
		ax.yaxis.set_major_locator(MaxNLocator(4))
		if chrom == legend_chrom:
			ax.legend(fontsize=7, loc='upper right', frameon=False)
		if ci == n_chroms - 1:
			ax.set_xlabel('Gene index', fontsize=8)

	fig.tight_layout(pad=1.0)
	for fmt in format:
		fpath = args.output + '.' + fmt
		fig.savefig(fpath, dpi=150)
		logger.info('Saved {}'.format(fpath))
	plt.close(fig)

	# Overall retention stats
	all_ref_genes = set()
	for chrom in d_chrom_paths:
		all_ref_genes.update(d_chrom_paths[chrom])
	if not args.count_duplicates:
		block_stats = compute_block_retention(syn_data, all_ref_genes,
											  args.ref, args.qry,
											  include_set, exclude_set)
	else:
		block_stats = compute_block_retention_dup(syn_data, all_ref_genes,
												   args.ref, args.qry,
												   include_set, exclude_set)
	stats_path = args.output + '.stats.tsv'
	_save_block_stats(block_stats, stats_path)
	logger.info('Saved {}'.format(stats_path))


# ---------------------------------------------------------------------------
#  parsing
# ---------------------------------------------------------------------------

def parse_collinearity(collinearity, ref, qry, min_same_block=25):
	"""Return {qry_sp: Counter(ref_gene -> multiplicity)}."""
	d_syn = {}
	for sp in qry:
		d_syn[sp] = Counter()

	qry_set = set(qry)
	for rc in XCollinearity(collinearity):
		if rc.chr1 == rc.chr2 and rc.N < min_same_block:
			continue
		sp1, sp2 = rc.species
		if sp1 == ref and sp2 in qry_set:
			for g1, g2 in rc.pairs:
				d_syn[sp2][g1] += 1
		elif sp2 == ref and sp1 in qry_set:
			for g1, g2 in rc.pairs:
				d_syn[sp1][g1] += 1
	return d_syn


def parse_gff(gff, ref):
	"""Return ({chrom: [(start, gene), ...]}, {chrom: [gene_names_in_order]})."""
	d_genes = OrderedDict()
	d_paths = OrderedDict()

	for line in XGff(gff):
		if line.species != ref:
			continue
		chrom = line.chrom
		if chrom not in d_genes:
			d_genes[chrom] = []
		d_genes[chrom].append((line.start, line.gene))

	for chrom, items in d_genes.items():
		items.sort(key=lambda x: x[0])
		path = [gene for _, gene in items]
		d_paths[chrom] = path

	return d_genes, d_paths


# ---------------------------------------------------------------------------
#  sliding window
# ---------------------------------------------------------------------------

def sliding_retention(path, syn, window_size, window_step,
					   include_set=None, exclude_set=None, count_duplicates=False):
	"""Sliding window retention rate along a gene path.

	syn: set (dedup mode) or Counter (count-duplicates mode).
	"""
	points = []
	for i in range(0, len(path), window_step):
		start, end = i, i + window_size
		if end > len(path):
			end = len(path)
		if end - start < window_size // 2:
			continue
		bin_genes = path[start:end]

		if include_set is not None:
			denom_genes = [g for g in bin_genes if g in include_set]
		elif exclude_set is not None:
			denom_genes = [g for g in bin_genes if g not in exclude_set]
		else:
			denom_genes = bin_genes

		n_denom = len(denom_genes)
		if n_denom == 0:
			continue

		if count_duplicates:
			n_syn = sum(syn.get(g, 0) for g in denom_genes)
		else:
			n_syn = sum(1 for g in denom_genes if g in syn)
		rate = n_syn / n_denom
		pos = (start + end) // 2
		points.append((pos, rate))
	return points


# ---------------------------------------------------------------------------
#  overall stats
# ---------------------------------------------------------------------------

def compute_block_retention(syn_sets, all_ref_genes, ref, qry,
							 include_set=None, exclude_set=None):
	"""Overall retention: |syntenic| / |all| (dedup mode)."""
	stats = []
	for qry_sp in qry:
		syn_g = syn_sets.get(qry_sp, set())
		if include_set is not None:
			all_g = all_ref_genes & include_set
		elif exclude_set is not None:
			all_g = all_ref_genes - exclude_set
		else:
			all_g = all_ref_genes
		n_all = len(all_g)
		n_syn = len(syn_g & all_g)
		rate = n_syn / n_all if n_all > 0 else 0
		stats.append((ref, qry_sp, n_syn, n_all, rate))
	return stats


def compute_block_retention_dup(syn_counters, all_ref_genes, ref, qry,
								 include_set=None, exclude_set=None):
	"""Overall retention with multiplicity (count-duplicates mode)."""
	stats = []
	for qry_sp in qry:
		counter = syn_counters.get(qry_sp, Counter())
		if include_set is not None:
			denom = all_ref_genes & include_set
		elif exclude_set is not None:
			denom = all_ref_genes - exclude_set
		else:
			denom = all_ref_genes
		n_all = len(denom)
		n_syn = sum(counter.get(g, 0) for g in denom)
		rate = n_syn / n_all if n_all > 0 else 0
		stats.append((ref, qry_sp, n_syn, n_all, rate))
	return stats


def _save_block_stats(stats, path):
	with open(path, 'w') as fout:
		fout.write('#ref\tqry\tsyntenic_genes\tall_genes\tretention\n')
		for ref_sp, qry_sp, n_syn, n_all, rate in stats:
			fout.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(
				ref_sp, qry_sp, n_syn, n_all, rate))


# ---------------------------------------------------------------------------
#  helpers
# ---------------------------------------------------------------------------

def _load_gene_list(path):
	"""Load gene IDs: split by whitespace, then by comma (for embedded lists)."""
	genes = set()
	with open(path) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			for token in line.split():
				for gene in token.split(','):
					gene = gene.strip()
					if gene:
						genes.add(gene)
	return genes
