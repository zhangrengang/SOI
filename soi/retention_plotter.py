# coding: utf-8
"""retention: gene retention rate along chromosomes.

Computes the proportion of reference genes that have syntenic orthologs
in each query species, using sliding-window and per-chromosome step plots
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

	syn_genes = parse_collinearity(args.collinearity, args.ref, args.qry,
									min_same_block=getattr(args, 'min_same_block', 25))
	if not args.count_duplicates:
		syn_sets = {}
		for qry, counter in syn_genes.items():
			syn_sets[qry] = set(counter.keys())
		syn_data = syn_sets
	else:
		syn_data = syn_genes

	d_chrom_genes, d_chrom_paths = parse_gff(args.gff, args.ref)
	include_set = _load_gene_list(args.include) if args.include else None
	exclude_set = _load_gene_list(args.exclude) if args.exclude else None

	if args.chrs:
		plot_chroms = [c for c in args.chrs if c in d_chrom_paths]
		if not plot_chroms:
			logger.error('None of the specified chromosomes found in GFF')
			sys.exit(1)
	else:
		candidates = [c for c in d_chrom_paths if len(d_chrom_paths[c]) >= args.min_genes]
		plot_chroms = sort_version(candidates)
	logger.info('Chromosomes to plot: {} ({} total, min_genes={})'.format(
		len(plot_chroms), len(d_chrom_paths), args.min_genes))

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

	# Overall retention per query
	all_ref_genes = set()
	for chrom in d_chrom_paths.values():
		all_ref_genes.update(chrom)
	if include_set is not None:
		denom_all = all_ref_genes & include_set
	elif exclude_set is not None:
		denom_all = all_ref_genes - exclude_set
	else:
		denom_all = all_ref_genes
	n_denom_all = len(denom_all)

	overall_rates = {}
	for qry in args.qry:
		if not args.count_duplicates:
			syn_g = syn_data.get(qry, set())
			n_syn = len(syn_g & denom_all)
		else:
			counter = syn_data.get(qry, Counter())
			n_syn = sum(counter.get(g, 0) for g in denom_all)
		overall_rates[qry] = n_syn / n_denom_all if n_denom_all > 0 else 0

	# Plot
	n_chroms = len(plot_chroms)
	figsize = (8, 0.8 * n_chroms)
	fig, axes = plt.subplots(n_chroms, 1, figsize=figsize,
							  sharex=True, sharey=True, squeeze=False)

	ymax = 0
	chrom_xmax = {}
	for qry in args.qry:
		for chrom in plot_chroms:
			pts = all_data.get(qry, {}).get(chrom, [])
			if pts:
				ymax = max(ymax, max(r for _, r in pts))
				chrom_xmax[chrom] = max(chrom_xmax.get(chrom, 0),
										max(p for p, _ in pts))
	if not args.count_duplicates:
		ymax = min(1.0, ymax * 1.05) if ymax > 0 else 1.0
	else:
		ymax = ymax * 1.05 if ymax > 0 else 1.0

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
			label = '{} ({:.2f})'.format(qry, overall_rates.get(qry, 0))
			ax.step(positions, rates, color=color, linewidth=0.8,
					label=label, where='mid')
		ax.set_ylim(0, ymax)
		ax.minorticks_on()
		ax.text(1.01, 0.5, chrom, transform=ax.transAxes, va='center',
				fontsize=9, rotation=90)
		ax.set_ylabel('Retention', fontsize=9)
		ax.yaxis.set_major_locator(MaxNLocator(4))
		ax.tick_params(labelsize=8)
		if ci == n_chroms - 1:
			ax.set_xlabel('Gene index', fontsize=9)

	# Legend
	fig.tight_layout(pad=0.8, h_pad=0.3)
	if len(args.qry) <= 3:
		_ax = axes[plot_chroms.index(legend_chrom)][0]
		_ax.legend(fontsize=8, loc='upper right', frameon=False)
	else:
		# Reserve right-side space and place legend in figure coordinates
		fig.subplots_adjust(right=0.78)
		h, l = [], []
		for ax in axes[:, 0]:
			_h, _l = ax.get_legend_handles_labels()
			if _h and not h:
				h, l = _h, _l
				break
		if h:
			fig.legend(h, l, fontsize=7, loc='center left',
					   frameon=False, bbox_to_anchor=(0.98, 0.5),
					   bbox_transform=fig.transFigure)

	for fmt in format:
		fpath = args.output + '.' + fmt
		fig.savefig(fpath, dpi=150)
		logger.info('Saved {}'.format(fpath))
	plt.close(fig)

	# Stats
	stats_path = args.output + '.stats.tsv'
	with open(stats_path, 'w') as fout:
		fout.write('#ref\tqry\tsyntenic_genes\tall_genes\tretention\n')
		for qry in args.qry:
			fout.write('{}\t{}\t.\t.\t{:.4f}\n'.format(
				args.ref, qry, overall_rates[qry]))
	logger.info('Saved {}'.format(stats_path))


# ---------------------------------------------------------------------------
#  parsing
# ---------------------------------------------------------------------------

def parse_collinearity(collinearity, ref, qry, min_same_block=25):
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
		d_paths[chrom] = [gene for _, gene in items]
	return d_genes, d_paths


def sliding_retention(path, syn, window_size, window_step,
					   include_set=None, exclude_set=None, count_duplicates=False):
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


def _load_gene_list(path):
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
