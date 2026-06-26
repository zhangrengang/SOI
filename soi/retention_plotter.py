# coding: utf-8
"""retention: gene retention rate along chromosomes.

Similar to depth (synteny depth / ploidy), but computes the proportion of
reference genes that have syntenic orthologs in each query species, using
sliding-window and per-chromosome line plots.
"""
import sys
import argparse
from math import ceil
from collections import OrderedDict
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from .mcscan import XCollinearity, XGff
from .RunCmdsMP import logger

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42


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
	# Normalize defaults
	if args.output is None:
		sps = [args.ref] + args.qry
		short = [x[:2] for x in sps]
		args.output = 'ret.' + '-'.join(short) + '_w' + str(args.window_size)
	format = args.format if isinstance(args.format, list) else [args.format]
	outfigs = [args.output + '.' + fmt for fmt in format]
	logger.info('{} x {} species'.format(len(args.qry), args.ref))

	# Load synteny -> syntenic gene pairs per query
	syn_genes = parse_collinearity(args.collinearity, args.ref, args.qry,
									min_same_block=getattr(args, 'min_same_block', 25))

	# Load GFF -> per-chromosome gene ordering for reference
	d_chrom_genes, d_chrom_paths = parse_gff(args.gff, args.ref)

	# Optional gene filtering
	include_set = _load_gene_list(args.include) if args.include else None
	exclude_set = _load_gene_list(args.exclude) if args.exclude else None

	# Determine chromosomes
	all_chroms = sorted(d_chrom_paths.keys())
	if args.chrs:
		plot_chroms = [c for c in args.chrs if c in d_chrom_paths]
		if not plot_chroms:
			logger.error('None of the specified chromosomes found in GFF')
			sys.exit(1)
	else:
		plot_chroms = [c for c in all_chroms if len(d_chrom_paths[c]) >= args.min_genes]

	logger.info('Chromosomes to plot: {} ({} total, min_genes={})'.format(
		len(plot_chroms), len(all_chroms), args.min_genes))

	# Sliding-window retention per query per chromosome
	all_data = OrderedDict()  # qry -> {chrom -> [(pos, rate), ...]}
	for qry in args.qry:
		syn_set = syn_genes.get(qry, set())
		chrom_data = OrderedDict()
		for chrom in plot_chroms:
			path = d_chrom_paths[chrom]
			points = sliding_retention(path, syn_set, args.window_size,
										args.window_step, include_set, exclude_set)
			if points:
				chrom_data[chrom] = points
		all_data[qry] = chrom_data

	# Plot
	n_species = len(args.qry)
	n_chroms = len(plot_chroms)
	ncols = n_species
	nrows = n_chroms

	figsize = (4 * ncols, 2.5 * nrows)
	fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
							  sharex='col', sharey='row', squeeze=False)

	# Global y-max for shared scale
	ymax = 0
	for qry in args.qry:
		for chrom in plot_chroms:
			pts = all_data.get(qry, {}).get(chrom, [])
			if pts:
				ymax = max(ymax, max(r for _, r in pts))
	ymax = min(1.0, ymax * 1.05) if ymax > 0 else 1.0

	for ci, chrom in enumerate(plot_chroms):
		for qi, qry in enumerate(args.qry):
			ax = axes[ci][qi]
			pts = all_data.get(qry, {}).get(chrom, [])
			if pts:
				positions = [p for p, _ in pts]
				rates = [r for _, r in pts]
				ax.plot(positions, rates, color='#2166ac', linewidth=0.8)
				ax.axhline(y=np.mean(rates), color='#b2182b', linewidth=0.5, linestyle='--')
			ax.set_ylim(0, ymax)
			ax.yaxis.set_major_locator(MaxNLocator(4))
			if ci == 0:
				ax.set_title(qry, fontsize=9)
			if qi == 0:
				ax.set_ylabel('{}\nRetention'.format(chrom), fontsize=8)
			if ci == nrows - 1:
				ax.set_xlabel('Gene index', fontsize=8)

	fig.tight_layout(pad=1.0)
	for fmt in format:
		fpath = args.output + '.' + fmt
		fig.savefig(fpath, dpi=150)
		logger.info('Saved {}'.format(fpath))
	plt.close(fig)

	# Per-block retention stats (denominator = all ref genes from GFF)
	all_ref_genes = set()
	for chrom in d_chrom_paths:
		all_ref_genes.update(d_chrom_paths[chrom])
	block_stats = compute_block_retention(syn_genes, all_ref_genes, args.ref, args.qry,
										  include_set, exclude_set)
	stats_path = args.output + '.stats.tsv'
	_save_block_stats(block_stats, stats_path)
	logger.info('Saved {}'.format(stats_path))


def parse_collinearity(collinearity, ref, qry, min_same_block=25):
	"""Return {qry_sp: set(ref_genes_with_synteny)}."""
	d_syn = {}
	for sp in qry:
		d_syn[sp] = set()

	qry_set = set(qry)
	for rc in XCollinearity(collinearity):
		if rc.chr1 == rc.chr2 and rc.N < min_same_block:
			continue
		sp1, sp2 = rc.species
		if sp1 == ref and sp2 in qry_set:
			for g1, g2 in rc.pairs:
				d_syn[sp2].add(g1)
		elif sp2 == ref and sp1 in qry_set:
			for g1, g2 in rc.pairs:
				d_syn[sp1].add(g1)
	return d_syn


def parse_gff(gff, ref):
	"""Return ({chrom: set(gene_names)}, {chrom: [gene_names_in_order]})."""
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


def sliding_retention(path, syn_set, window_size, window_step,
					   include_set=None, exclude_set=None):
	"""Sliding window retention rate along a gene path.

	Retention = |syntenic_genes_in_window| / |total_genes_in_window|.
	Returns [(gene_index, rate), ...].
	"""
	points = []
	for i in range(0, len(path), window_step):
		start, end = i, i + window_size
		if end > len(path):
			end = len(path)
		if end - start < window_size // 2:
			continue
		bin_genes = path[start:end]

		# Filter denominator
		if include_set is not None:
			denom_genes = [g for g in bin_genes if g in include_set]
		elif exclude_set is not None:
			denom_genes = [g for g in bin_genes if g not in exclude_set]
		else:
			denom_genes = bin_genes

		n_denom = len(denom_genes)
		if n_denom == 0:
			continue

		n_syn = sum(1 for g in denom_genes if g in syn_set)
		rate = n_syn / n_denom
		pos = (start + end) // 2
		points.append((pos, rate))
	return points


def compute_block_retention(syn_genes, all_ref_genes, ref, qry, include_set=None, exclude_set=None):
	"""Compute overall retention: syntenic ref genes / all ref genes."""
	stats = []
	for qry_sp in qry:
		syn_g = syn_genes.get(qry_sp, set())
		if include_set is not None:
			all_g = all_ref_genes & include_set
			syn_g = syn_g & include_set
		elif exclude_set is not None:
			all_g = all_ref_genes - exclude_set
			syn_g = syn_g - exclude_set
		else:
			all_g = all_ref_genes
		n_all = len(all_g)
		n_syn = len(syn_g & all_g)  # only count syntenic genes that are also in the denominator
		rate = n_syn / n_all if n_all > 0 else 0
		stats.append((ref, qry_sp, n_syn, n_all, rate))
	return stats


def _save_block_stats(stats, path):
	"""Save per-block retention to TSV."""
	with open(path, 'w') as fout:
		fout.write('#ref\tqry\tsyntenic_genes\tall_genes\tretention\n')
		for ref_sp, qry_sp, n_syn, n_all, rate in stats:
			fout.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(
				ref_sp, qry_sp, n_syn, n_all, rate))


def _load_gene_list(path):
	"""Load gene IDs from file, one per line."""
	genes = set()
	with open(path) as f:
		for line in f:
			line = line.strip()
			if line and not line.startswith('#'):
				genes.add(line)
	return genes
