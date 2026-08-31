import sys
import argparse
import numpy as np
import collections
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from lazy_property import LazyWritableProperty as lazyproperty
from .colors import Colors
from .mcscan import Collinearity, Gff, XCollinearity
from .RunCmdsMP import logger
from .style import style_species


def evaluate_args(parser):
	parser.add_argument('-s', '-synteny', required=True, nargs='+',
						dest='collinearities', metavar='FILE',
						help='Collinearity files [required]')
	parser.add_argument('-o', '-orthology', nargs='+', default=None,
						dest='orthologs', metavar='FILE',
						help='Ortholog files for OI [optional; OI panels hidden if omitted]')
	parser.add_argument('-g', '-gff', required=True, nargs='+',
						dest='gff', metavar='GFF',
						help='GFF file [required]')
	parser.add_argument('-r', '-ref', required=True, type=str,
						dest='ref', metavar='SPECIES',
						help='Reference species for fractionation rate calculation [required]')
	parser.add_argument('-q', '-qry', type=str, nargs='+', default=None,
						dest='qry', metavar='SPECIES',
						help='Target query species to evaluate (only pairs involving these species)')
	parser.add_argument('-pre', '-prefix', type=str, default=None,
						dest='pre', metavar='PREFIX',
						help='Output prefix [default: <synteny>_<ref>]')
	parser.add_argument('--blocks', metavar='NAME', type=str, nargs='*', default=None,
						dest='blocks',
						help='treat -s input(s) as paralog blocks file(s) and aggregate by branch; optional branch names, one-to-one with -s files (default: parsed from file names <prefix>.<branch>.blocks) [default: off]')
	parser.add_argument('--figsize', type=float, nargs=2, default=(10, 12),
						dest='figsize', metavar=('W', 'H'),
						help='Figure size in inches [default: 10 12]')


def evaluate_main(**kargs):
	import os
	if kargs.get('pre') is None:
		basename = os.path.splitext(os.path.basename(kargs['collinearities'][0]))[0]
		kargs['pre'] = '{}_{}'.format(basename, kargs['ref'])
	eval(**kargs)


def eval(collinearities, orthologs, gff, qry=None, ref=None, pre=None, figsize=(10, 12),
		 blocks=None):
	# --blocks: aggregate by branch; names one-to-one with -s files or parsed from file names
	d_file_branch = {}
	if blocks is not None:
		from .ks_plotter import parse_branch
		files = [collinearities] if isinstance(collinearities, str) else collinearities
		if blocks:
			if len(blocks) != len(files):
				logger.error('--blocks names ({}) must match -s files ({}) one-to-one'.format(
					len(blocks), len(files)))
				sys.exit(1)
			for f, b in zip(files, blocks):
				d_file_branch[f] = b
		else:  # no names: parse from file names
			for f in files:
				d_file_branch[f] = parse_branch(f)
	d_rcs = {}
	d_refgenes = {}
	for rc in XCollinearity(collinearities, orthologs=orthologs, gff=gff):
		rc.fr = rc.fractionation_rate(ref=ref)
		if rc.fr is None:
			continue
		sp1, sp2 = rc.species1, rc.species2
		genes1, genes2 = rc.genes
		branch = None
		if blocks is not None:
			branch = d_file_branch.get(getattr(rc, 'source_file', None))
		elif sp1 == sp2:
			continue
		if qry and sp1 not in qry and sp2 not in qry:
			continue
		if sp2 == ref:
			sp1, sp2 = sp2, sp1
			genes1, genes2 = genes2, genes1
		xrc = SynRec(rc)
		key = (branch, sp1, sp2) if blocks is not None else (sp1, sp2)
		try:
			d_rcs[key] += [xrc]
		except KeyError:
			d_rcs[key] = [xrc]
		try:
			d_refgenes[key] += genes1
		except KeyError:
			d_refgenes[key] = genes1

	for spp, rcs in d_rcs.items():
		rcs = SynRecs(rcs)
		counts = collections.Counter(d_refgenes[spp])
		rcs.refcounts = np.array(
			sorted(collections.Counter(counts.values()).items()))
		d_rcs[spp] = rcs
	pre = pre or 'evaluate'
	outfig = pre  # without extension, plot_eval appends .png and .pdf
	plot_eval(d_rcs, outfig, figsize=figsize)
	logger.info('Saved: %s.png, %s.pdf', outfig, outfig)



def plot_eval(d_rcs, outfig, legend_fontsize=9, figsize=(10, 12)):
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(4, 7)

    ax0 = fig.add_subplot(gs[:, 6])  # legend
    ax1 = fig.add_subplot(gs[0, 0:3])  # accumulation N
    ax2 = fig.add_subplot(gs[1, 0:3])  # decay N
    ax3 = fig.add_subplot(gs[2, 0:3])  # histo of fr
    ax4 = fig.add_subplot(gs[0, 3:6])  # dots
    ax5 = fig.add_subplot(gs[1, 3:6])  # dots of 50
    ax6 = fig.add_subplot(gs[2, 3:6])  # bar
    ax7 = fig.add_subplot(gs[3, 0:3])
    ax8 = fig.add_subplot(gs[3, 3:6])

    n = len(d_rcs)
    colors = Colors(n).colors
    has_oi = all(rcs.has_oi for rcs in d_rcs.values())
    line = ['species1', 'species2',  'block number', 'gene number', 'min size',
            'max size', 'median size', 'mean size', 'size 50',
            'FR50', 'OI50']
    logger.info('\t'.join(line))
    for i, (key, rcs) in enumerate(sorted(d_rcs.items(), key=str)):
        if len(key) == 3:  # blocks mode: (branch, sp1, sp2)
            from .style import format_species_label
            branch, sp1, sp2 = key
            s1, _ = format_species_label(sp1)
            s2, _ = format_species_label(sp2)
            _sp = '{}: {} {}'.format(branch, s1, s2)
        else:
            branch = None
            sp1, sp2 = key
            _sp = convert_sp(sp2)
        i50, sn50, fn50 = rcs.xn50  # index 50, block size 50, frac rate 50
        fm, sm = np.median(rcs.extended_frs), np.median(rcs.extended_ns)
        # frac rate 50, block size 50
        om = np.median(rcs.extended_ois) if has_oi else None  # OI 50
        ns = rcs.ns  # block sizes: Ns
        nb = len(ns)  # block number
        bm = sum(rcs.refcounts[:, 1])  # gene number of all blocks
        kargs = dict(color=colors[i], alpha=1, label=_sp)
        ax0.plot(-2, -2, linestyle='-', marker='o',  **kargs)

        ax1.plot(range(1, nb+1), np.cumsum(ns),
                 drawstyle="steps-post",  **kargs)
        ax2.plot(range(1, nb+1), ns, drawstyle="steps-post", **kargs)
        ax3.scatter(len(rcs.ns), bm, **kargs)
        hist_plot(rcs.extended_frs, ax4, bins=40, **kargs)
        ax5.scatter(fm, sm, **kargs)

        if has_oi:
            ax6.scatter(fm, om, **kargs)
            hist_plot(rcs.extended_ois, ax7, bins=40, **kargs)
        ax8.plot(rcs.refcounts[:, 0], rcs.refcounts[:, 1], **kargs)

        line = [branch if branch is not None else sp1, sp2, nb, bm,
                int(min(ns)), int(max(ns)), int(np.median(ns)),
                round(np.mean(ns), 1), sm, round(fm, 2),
                round(om, 2) if has_oi else 'NA']
        logger.info('\t'.join(map(str, line)))
    # ax
    # 1	4
    # 2	5
    # 3	6
    # 7	8
    for ax, xlab, ylab in (
            [ax1, '# of blocks', 'Cumulative number of genes'],
            [ax2, '# of blocks', 'Block size'],
            [ax3, 'Total number of blocks', 'Total number of genes'],
            [ax4, 'Fractionation rate', 'Number of genes'],
            [ax5, 'Fractionation rate 50', 'Block size 50'],
            [ax6, 'Fractionation rate 50', 'Orthology index 50'],
            [ax7, 'Orthology index', 'Number of genes'],
            [ax8, 'Copy number', 'Number of genes'],
    ):
        if ax in (ax6, ax7) and not has_oi:
            ax.axis('off')  # OI panels hidden without orthologs input
            continue
        set_labels(ax, xlab, ylab)

    # legend
    ncols = n // 40 + 1
    ax0.set_xlim(0, 1)
    ax0.set_ylim(0, 1)
    [style_species(t) for t in ax0.legend(loc=(-1.5, 0), fancybox=False,
                                          frameon=False, ncols=ncols).get_texts()]
    ax0.xaxis.set_tick_params(length=0)
    ax0.spines['right'].set_color('none')
    ax0.spines['top'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.xaxis.set_ticks([])
    ax0.yaxis.set_ticks([])
# ax0.set_title('Species', loc='left')

    plt.subplots_adjust(hspace=0.3, wspace=1.8)
    for fmt in ('png', 'pdf'):
        plt.savefig('{}.{}'.format(outfig, fmt), bbox_inches='tight')


def hist_plot(data, ax, bins=10, alpha=1, **kargs):
    n, bins, patches = ax.hist(data, bins=bins, alpha=0, **kargs)
    Xs, Ys = [], []
    for i in range(len(bins)-1):
        X = (bins[i] + bins[i+1])/2
        Xs.append((bins[i] + bins[i+1])/2)
        Ys.append(n[i])

    ax.plot(Xs, Ys, alpha=alpha, **kargs)


def set_labels(ax, xlab, ylab, **kargs):
    ax.set_xlabel(xlab, **kargs)
    ax.set_ylabel(ylab, **kargs)


def set_legends(axs, **kargs):
    for ax in axs:
        ax.legend(**kargs)


def convert_sp(sp):
    return sp


class SynRec:
    def __init__(self, rc):
        self.N = rc.N
        self.fr = rc.fr
        self.oi = getattr(rc, 'oi', None)  # None if no orthologs input


class SynRecs:
    def __init__(self, rcs):
        self.rcs = rcs

    @property
    def has_oi(self):
        return any(rc.oi is not None for rc in self.rcs)

    @lazyproperty
    def matrix(self):
        data = []
        for rc in self.rcs:
            oi = rc.oi if rc.oi is not None else 0  # placeholder; panels hidden without OI
            data += [[rc.N, rc.fr, oi]]
        return np.array(sorted(data, key=lambda x: -x[0]))

    @lazyproperty
    def ns(self):
        return self.matrix[:, 0]

    @lazyproperty
    def frs(self):
        return self.matrix[:, 1]

    @lazyproperty
    def ois(self):
        return self.matrix[:, 2]

    def get_nx0(self, cutoff=50):
        self.size = sum(self.ns)
        accum = 0
        for i, x in enumerate(self.ns):
            accum += x
            if 1e2 * accum / self.size >= cutoff:
                return x, i+1

    @lazyproperty
    def xn50(self):
        sn50, i50 = self.get_nx0()
        fn50 = list(sorted(self.frs, reverse=1))[i50 - 1]
        return i50, sn50, fn50

    @lazyproperty
    def extended_frs(self):
        frs = []
        for n, fr, oi in self.matrix:
            frs += [fr] * int(n)
        return frs

    @lazyproperty
    def extended_ns(self):
        nrs = []
        for n in self.ns:
            nrs += [n] * int(n)
        return nrs

    @lazyproperty
    def extended_ois(self):
        ois = []
        for n, fr, oi in self.matrix:
            ois += [oi] * int(n)
        return ois
