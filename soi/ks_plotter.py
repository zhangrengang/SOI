import sys,os
import argparse
from collections import OrderedDict
import numpy as np
import matplotlib as mpl

from .style import format_species_title
from .colors import Colors
from .RunCmdsMP import logger

def add_ksplot_args(parser):
	parser.add_argument('--kaks', metavar='FILE', type=str, nargs='+', default=None,
		help='KaKs input file(s). Formats: WGDI, or KaKsCalculator/yn00/4DTV whose '
			 'Sequence IDs are labelled SPECIES|GENE-SPECIES|GENE (curves are '
			 'grouped by the species pair parsed from those IDs) [required]')
	parser.add_argument('--blocks', metavar='FILE', type=str, nargs='+', default=None,
		help='Paralog blocks file(s), Ks joined from --kaks [optional]')
	parser.add_argument('--block-branches', metavar='NAME', type=str, nargs='+', default=None,
		dest='block_branches', help='Branch name(s), one-to-one with --blocks [default: parsed from file names]')
	parser.add_argument('--pair', metavar='FILE', type=str, default=None,
		help='Species pair control file [optional]')
	parser.add_argument('-pre', '--prefix', metavar='PREFIX', type=str, default=None,
		dest='figure', help='Output prefix [default: --pair or --kaks + suffix + .pdf]')
	parser.add_argument('--max-ks', metavar='FLOAT', type=float, default=2,
		help='Max Ks [default=%(default)s]')
	parser.add_argument('--bins-per-ks', metavar='INT', type=int, default=50,
		help='Bins per Ks unit [default=%(default)s]')
	parser.add_argument('--normed', action='store_true', default=False,
		help='Normalize histogram [default=%(default)s]')
	parser.add_argument('--yn00', action='store_true', default=False,
		help='Use YN00 method [default=%(default)s]')
	parser.add_argument('--method', metavar='STR', type=str, default='NG86',
		help='KaKs method [default=%(default)s]')
	parser.add_argument('--fdtv', action='store_true', default=False,
		help='Use 4DTV [default=%(default)s]')
	parser.add_argument('--xlabel', metavar='STR', type=str, default='Ks',
		help='X-axis label [default=%(default)s]')
	parser.add_argument('--homology-class', metavar='FILE', type=str, default=None,
		help='Homology class file [optional]')
	parser.add_argument('--only-class', action='store_true', default=False,
		help='Only plot homology in --homology-class [default=%(default)s]')
	parser.add_argument('--plot-others', action='store_true', default=False,
		help='Plot homology not in --homology-class as a single line [default=%(default)s]')
	parser.add_argument('--out-data', action='store_true', default=False,
		help='Output plot data [default=%(default)s]')
	parser.add_argument('--remove-same-chr', action='store_true', default=False,
		help='Remove data on the same chromosome [default=%(default)s]')
	parser.add_argument('--gff', type=str, default=None,
		help='MCScanX GFF file for --remove-same-chr [default=%(default)s]')
	parser.add_argument('--plot-type', '-p', metavar='TYPE', type=str,
		nargs='+', default=['all'],
		choices=['hist', 'density', 'ridge', 'all'],
		help='Plot type(s): hist, density, ridge, all. [default=%(default)s]')


def ksplot_args(parser):
	add_ksplot_args(parser)


def xmain(**kargs):
	main(**kargs)


def main(kaks=None, pair=None, blocks=None, block_branches=None, figure=None,
		 max_ks=2, bins_per_ks=50, normed=False, yn00=False, method='NG86',
		 fdtv=False, xlabel='Ks', homology_class=None, only_class=False,
		 plot_others=False, out_data=False, remove_same_chr=False,
		 gff=None, plot_type=None, **kargs):

	if plot_type is None:
		plot_type = ['all']
	if 'all' in plot_type:
		plot_type = ['hist', 'density', 'ridge']

	# Auto figure name
	if fdtv and xlabel == 'Ks':
		xlabel = '4dTV'
	suffix = xlabel + '.yn00' if yn00 else xlabel
	if isinstance(kaks, str):
		kaks = [kaks]
	if blocks is not None:
		if isinstance(blocks, str):
			blocks = [blocks]
		if figure is None:
			if len(blocks) > 1:
				logger.warning('multiple --blocks files without -pre, using the first file name as prefix') # [hidden]
			figure = blocks[0].rsplit('.blocks', 1)[0] + '.' + suffix
	elif figure is None:
		if pair is not None:
			figure = pair + '.' + suffix
		else:
			if len(kaks) > 1:
				logger.warning('multiple --kaks files without -pre, using the first file name as prefix') # [hidden]
			figure = kaks[0] + '.' + suffix

	pairs = parse_pair(pair)
	homology_class_map = parse_homology_class(homology_class)
	d_gff = parse_gff(gff) if remove_same_chr else None

	if blocks is not None:
		if kaks is None:
			logger.error('--blocks requires --kaks (Ks source)')
			sys.exit(1)
		if homology_class_map is not None:
			logger.error('--blocks and --homology-class are mutually exclusive (branch is the class)')
			sys.exit(1)
		d_ks = parse_blocks_ks(blocks, kaks, block_branches=block_branches,
			pairs=pairs, max_ks=max_ks, fdtv=fdtv, yn00=yn00, method=method)
		pairs = None  # order keys are (sp1, sp2, branch); plot falls back to sorted keys
	else:
		d_ks = parse_kaks(kaks=kaks, pairs=pairs, max_ks=max_ks, fdtv=fdtv,
			yn00=yn00, method=method, d_gff=d_gff,
			homology_class=homology_class_map, only_class=only_class,
			plot_others=plot_others)

	if homology_class_map is not None and pairs is not None:
		_porder = {(p[0], p[1]): i for i, p in enumerate(pairs)}
		pairs = sorted(sorted(d_ks.keys()),
					   key=lambda x: _porder.get((x[0], x[1]), len(_porder)))

	if out_data:
		out_data_file = figure + '.data.xls'
		with open(out_data_file, 'w') as f:
			print('#parameters:', kargs, file=f)
	else:
		out_data_file = None

	nbins = int(max_ks * bins_per_ks)

	# one colour per species pair, shared by all plot types (keyed by pair,
	# not by plotting position: density/ridge drop groups with <3 points).
	# With a grouped --pair file (3rd column) pairs of the same group share
	# a colour; ungrouped pairs each keep their own.
	if pairs is not None:
		_order = list(pairs)
	else:
		_order = sorted(d_ks.keys())
	_gkeys = [p[2] if len(p) > 2 else p for p in _order]
	_groups = []
	for _g in _gkeys:
		if _g not in _groups:
			_groups.append(_g)
	_gcolor = {g: c for g, c in zip(_groups, Colors(len(_groups)).colors)}
	d_color = {tuple(p): _gcolor[g] for p, g in zip(_order, _gkeys)}

	# Draw selected plot types
	for ptype in plot_type:
		if ptype == 'hist':
			outfig = figure + '.hist.pdf'
			plot_hist(d_ks, outfig, order=pairs, d_color=d_color, normed=normed,
					  nbins=nbins, xlabel=xlabel, out_data=out_data_file,
					  max_ks=max_ks)
		elif ptype == 'density':
			outfig = figure + '.density.pdf'
			plot_density(d_ks, outfig, order=pairs, d_color=d_color, xlabel=xlabel,
						 max_ks=max_ks)
		elif ptype == 'ridge':
			outfig = figure + '.ridge.pdf'
			plot_ridge(d_ks, outfig, order=pairs, d_color=d_color, xlabel=xlabel,
					   max_ks=max_ks)


# =====================
# Plot functions
# =====================

def plot_hist(d_data, outFig, order=None, normed=False, xlabel='Ks',
			  nbins=300, ylabel=' of syntenic gene pairs', out_data=None,
			  max_ks=3, d_color=None):
	"""Histogram-based line plot (original behavior)."""
	import matplotlib.pyplot as plt

	if order is None:
		order = sorted(d_data.keys())
	if normed:
		ylabel = 'Percent' + ylabel
	else:
		ylabel = 'Number' + ylabel

	fig, ax = plt.subplots()
	x_list = []
	y_list = []
	pair_list = []
	pair_keys = []
	label_styles = []

	print('\t'.join(['species1', 'species2', 'category', 'genePairNumber',
					 'mean', 'median', 'peak', '95% CI']))
	for pair in order:
		data = np.array(d_data[pair])
		label, lstyle = generate_label(pair)
		pair_list.append(label)
		pair_keys.append(tuple(pair))
		label_styles.append(lstyle)
		data2 = np.concatenate([[0, max_ks], data])
		n, bins, patches = ax.hist(data2, bins=nbins, density=normed,
								   label=label)
		y_hist = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]
		x_list.append(y_hist)
		y_list.append(n)
		y_hist_fmt = format_ks(y_hist, nbins)
		peak_ks = max(zip(y_hist_fmt, n), key=lambda x: x[1])[0]
		if out_data:
			with open(out_data, 'a') as f:
				print('\t'.join(list(pair) + ['x'] + list(map(str, y_hist_fmt))), file=f)
				print('\t'.join(list(pair) + ['y'] + list(map(str, n))), file=f)
		tile2_5 = np.percentile(data, 2.5)
		tile97_5 = np.percentile(data, 97.5)
		ci95 = '{:.3f}-{:.3f}'.format(abs(tile2_5), tile97_5)
		category = '-' if len(pair) == 2 else pair[2]
		print('\t'.join(map(str, [pair[0], pair[1], category, len(data),
								  np.mean(data), np.median(data), peak_ks, ci95])))

	ax.cla()
	for i in range(len(pair_list)):
		ax.plot(x_list[i], y_list[i], label=pair_list[i], linewidth=1.5,
				color=(d_color or {}).get(pair_keys[i]))
	for t, s in zip(ax.legend(fontsize='x-small').get_texts(), label_styles):
		t.set_fontstyle(s)
	ax.minorticks_on()
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if str(outFig).endswith('.pdf'):
		plt.savefig(str(outFig)[:-4] + '.png', dpi=300)
	plt.savefig(outFig)
	plt.close()


def plot_density(d_data, outFig, order=None, xlabel='Ks', max_ks=3, d_color=None):
	"""Density plot using scipy.stats.gaussian_kde (pure Python, no R)."""
	from scipy.stats import gaussian_kde
	import matplotlib.pyplot as plt

	if order is None:
		order = sorted(d_data.keys())

	fig, ax = plt.subplots()
	xs = np.linspace(0, max_ks, 500)

	label_styles = []
	for pair in order:
		data = np.array(d_data[pair])
		data = data[(data > 0) & (data <= max_ks)]
		if len(data) < 3:
			continue
		try:
			kde = gaussian_kde(data)
			ys = kde(xs)
		except np.linalg.LinAlgError:
			continue
		label, lstyle = generate_label(pair)
		label_styles.append(lstyle)
		# Distinguish orthologs vs paralogs by line style
		if len(set(pair[:2])) == 1:
			ax.plot(xs, ys, label=label, linewidth=1.5, linestyle='--',
					color=(d_color or {}).get(tuple(pair)))
		else:
			ax.plot(xs, ys, label=label, linewidth=1.5,
					color=(d_color or {}).get(tuple(pair)))

	ax.set_xlabel(xlabel)
	ax.set_ylabel('Density')
	ax.set_xlim(0, max_ks)
	for t, s in zip(ax.legend(fontsize='x-small').get_texts(), label_styles):
		t.set_fontstyle(s)
	ax.minorticks_on()
	if str(outFig).endswith('.pdf'):
		plt.savefig(str(outFig)[:-4] + '.png', dpi=300)
	plt.savefig(outFig)
	plt.close()


def plot_ridge(d_data, outFig, order=None, xlabel='Ks', max_ks=3, d_color=None):
	"""Ridge plot (山岭图) — stacked density curves per species pair."""
	from scipy.stats import gaussian_kde
	import matplotlib.pyplot as plt

	if order is None:
		order = sorted(d_data.keys())

	# Precompute KDE for each group
	curves = []
	labels = []
	pair_keys = []
	xs = np.linspace(0, max_ks, 500)
	for pair in order:
		data = np.array(d_data[pair])
		data = data[(data > 0) & (data <= max_ks)]
		if len(data) < 3:
			continue
		try:
			kde = gaussian_kde(data)
			ys = kde(xs)
			ys = ys / ys.max()  # normalize to [0, 1] for visual comparison
		except np.linalg.LinAlgError:
			continue
		curves.append(ys)
		labels.append(generate_label(pair))  # (text, style)
		pair_keys.append(tuple(pair))

	if not curves:
		print("Warning: no valid data for ridge plot", file=sys.stderr)
		return

	n = len(curves)
	# row height shrinks with the number of groups so the figure stays
	# printable (~<=12 in tall); labels keep a fixed 8 pt, so their size
	# relative to the figure grows with n instead of shrinking when the
	# tall figure is scaled down for a page
	row_h = 1.0 if n <= 10 else max(0.35, 11.0 / n)
	fig, axes = plt.subplots(n, 1, sharex=True,
							 figsize=(8, row_h * n + 1),
							 gridspec_kw={'hspace': 0})
	if n == 1:
		axes = [axes]

	for i, (ax, ys, label) in enumerate(zip(axes, curves, labels)):
		text, lstyle = label
		color = (d_color or {}).get(pair_keys[i])
		ax.fill_between(xs, ys, alpha=0.5, color=color)
		ax.plot(xs, ys, color='black', linewidth=0.5)
		ax.set_yticks([])
		ax.set_ylim(0, 1.05)
		ax.set_xlim(0, max_ks)
		# Label on the right
		ax.text(1.01, 0.5, text, transform=ax.transAxes,
				va='center', ha='left', fontsize=8, fontstyle=lstyle)
		# Remove spines except bottom
		for spine in ['left', 'right', 'top']:
			ax.spines[spine].set_visible(False)

	axes[-1].set_xlabel(xlabel)
	axes[-1].spines['bottom'].set_visible(True)

	# Hide bottom spines for upper axes
	for ax in axes[:-1]:
		ax.spines['bottom'].set_visible(False)
		ax.tick_params(bottom=False)

	plt.tight_layout()
	plt.subplots_adjust(hspace=0)
	if str(outFig).endswith('.pdf'):
		plt.savefig(str(outFig)[:-4] + '.png', dpi=300)
	plt.savefig(outFig)
	plt.close()


# =====================
# Data parsing (unchanged)
# =====================

def parse_branch(block_file):
	'''branch name from blocks file name: <prefix>.<branch>.blocks -> <branch>'''
	base = os.path.basename(block_file)
	if not base.endswith('.blocks'):
		logger.error('--blocks file name must end with .blocks: {}'.format(block_file))
		sys.exit(1)
	base = base[:-len('.blocks')]
	if '.' in base:
		base = base.split('.', 1)[1]  # strip prefix
	return base


def parse_blocks_ks(blocks, kaks, block_branches=None, pairs=None, max_ks=None,
				   fdtv=False, yn00=False, method='NG86'):
	'''Join blocks pairs with Ks values; d_ks[(sp1, sp2, branch)] = [ks, ...]'''
	from .mcscan import XCollinearity, XKaKsParser

	# Ks lookup: bidirectional pair -> ks
	d_pair = {}
	for rc in XKaKsParser(kaks):
		rc.parse_pair()
		if rc.ks is None or rc.ks <= 0:
			continue
		g1, g2 = rc.pair
		d_pair[(g1, g2)] = rc.ks
		d_pair[(g2, g1)] = rc.ks
	logger.info('{} Ks pairs from {} kaks file(s)'.format(len(d_pair) // 2, kaks))

	if block_branches is not None:
		if len(block_branches) != len(blocks):
			logger.error('--block-branches ({}) must match --blocks ({}) one-to-one'.format(
				len(block_branches), len(blocks)))
			sys.exit(1)
		d_branch = dict(zip(blocks, block_branches))
	else:
		d_branch = {b: parse_branch(b) for b in blocks}

	if pairs is not None:
		pairs = set(pairs)

	d_ks = {}
	d_miss = {}
	for rc in XCollinearity(blocks):
		branch = d_branch.get(rc.source_file, rc.source_file)
		sp1, sp2 = rc.species
		for g1, g2 in rc.pairs:
			ks = d_pair.get((g1, g2))
			if ks is None:
				d_miss[branch] = d_miss.get(branch, 0) + 1
				continue
			if max_ks is not None and ks > max_ks:
				continue
			if pairs is not None and (sp1, sp2) not in pairs and (sp2, sp1) not in pairs:
				continue
			Xpair = (sp1, sp2, branch)
			try:
				d_ks[Xpair] += [ks]
			except KeyError:
				d_ks[Xpair] = [ks]
	for branch, n_miss in sorted(d_miss.items()):
		logger.warning('branch {}: {} pairs without Ks, skipped'.format(branch, n_miss))
	if not d_ks:
		logger.error('No Ks values joined from blocks; check --kaks and --blocks overlap')
		sys.exit(1)
	return d_ks


def parse_homology_class(homology_class):
	if homology_class is None:
		return None
	d = {}
	for line in open(homology_class):
		temp = line.strip().split()
		pair = tuple(temp[:2])
		try:
			CLASS = temp[2]
		except IndexError:
			CLASS = None
		if CLASS is not None:
			try:
				d[pair] += [CLASS]
			except KeyError:
				d[pair] = [CLASS]
		else:
			d[pair] = []
	return d


def parse_pair(pair):
	'''Read a species-pair control file: 'sp1 sp2 [group]' per line.

	A third column groups pairs (same group -> same colour in the plots
	and the group name is used as the curve label/category).
	'''
	if pair is None:
		return None
	pairs = []
	for line in open(pair):
		if line.startswith('#'):
			continue
		temp = line.strip().split()
		if len(temp) < 2:
			continue
		pairs.append(tuple(temp[:3]))  # (sp1, sp2[, group])
	return pairs


def parse_gff(gff):
	d = {}
	for line in open(gff):
		line = line.strip().split()
		chr, g = line[:2]
		d[g] = chr
	return d


def parse_kaks(kaks, pairs=None, max_ks=None, fdtv=False, yn00=False,
			   method='NG86', homology_class=None, only_class=False,
			   plot_others=False, d_gff=None):
	from .mcscan import KaKs, XKaKsParser

	if pairs is not None:
		# key by the first two columns so a control file with an optional
		# group column still matches; the full tuple (with group) is the key
		d_pairs = {}
		for _p in pairs:
			d_pairs[(_p[0], _p[1])] = _p
			d_pairs[(_p[1], _p[0])] = _p
	d_ks = {}
	for rc in XKaKsParser(kaks):
		rc.parse_pair()
		sp1, sp2 = rc.species
		g1, g2 = rc.pair
		if d_gff:
			if d_gff[g1] == d_gff[g2]:
				continue
		if pairs is None:
			pair = (sp1, sp2)
		else:
			pair = d_pairs.get((sp1, sp2))
			if pair is None:
				continue
		if rc.ks is None or rc.ks <= 0:
			continue
		if max_ks is not None and rc.ks > max_ks:
			continue

		Xpairs = [pair]
		if homology_class is not None:
			if only_class:
				Xpairs = []
			else:
				Xpairs = [tuple(list(pair) + ['ALL'])]
			if (g1, g2) in homology_class:
				classes = homology_class[(g1, g2)]
			elif (g2, g1) in homology_class:
				classes = homology_class[(g2, g1)]
			else:
				if only_class:
					continue
				elif plot_others:
					classes = ['OTHERS']
				else:
					classes = []
			for CLASS in classes:
				Xpairs.append(tuple(list(pair) + [CLASS]))
		for pair in Xpairs:
			try:
				d_ks[pair] += [rc.ks]
			except KeyError:
				d_ks[pair] = [rc.ks]
	return d_ks


# =====================
# Utilities (unchanged)
# =====================

def generate_label(pair):
	if len(pair) > 2:
		from .style import format_species_label
		text, lstyle = format_species_label(pair[2])
		return (text, lstyle)
	return format_species_title('-'.join(pair), sep='-')


def format_ks(values, bins):
	return [int(round(v * bins, 0)) * 1.0 / bins for v in values]


if __name__ == '__main__':
	main()
