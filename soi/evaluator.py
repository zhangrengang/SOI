import sys
import numpy as np
import collections
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42

from lazy_property import LazyWritableProperty as lazyproperty

from .mcscan import Collinearity, Gff, XCollinearity
from .colors import Colors

def eval(collinearities, orthologs, gff, ref, pre=None):
	d_rcs = {}
	d_refgenes = {}
	for rc in XCollinearity(collinearities, orthologs=orthologs, gff=gff):
		rc.fr = rc.fractionation_rate(ref=ref)
		if rc.fr is None:
			continue
		sp1, sp2 = rc.species1, rc.species2
		genes1, genes2 = rc.genes
		if sp1 == sp2:
			continue
		if sp2 == ref:
			sp1, sp2 = sp2, sp1
			genes1, genes2 = genes2, genes1
		xrc = SynRec(rc)
		try: d_rcs[sp2] += [xrc]
		except KeyError: d_rcs[sp2] = [xrc]
		try: d_refgenes[sp2] += genes1
		except KeyError: d_refgenes[sp2] = genes1

	for sp, rcs in d_rcs.items():
		rcs = SynRecs(rcs)
		counts = collections.Counter(d_refgenes[sp])
		rcs.refcounts = np.array(sorted(collections.Counter(counts.values()).items()))
		d_rcs[sp] = rcs
	outfig = 'test.png'
	plot_eval(d_rcs, outfig)
def main():
	collinearities, orthologs, gff, ref = sys.argv[1:5]
	eval(collinearities, orthologs, gff, ref)
def plot_eval(d_rcs, outfig, legend_fontsize=9):
	fig = plt.figure(figsize=(10,12))
	gs = GridSpec(4, 7)
	
	ax0 = fig.add_subplot(gs[:, 6]) # legend
	ax1 = fig.add_subplot(gs[0, 0:3]) # accumulation N
	ax2 = fig.add_subplot(gs[1, 0:3]) # decay N
	ax3 = fig.add_subplot(gs[2, 0:3]) # histo of fr
	ax4 = fig.add_subplot(gs[0, 3:6]) # dots
	ax5 = fig.add_subplot(gs[1, 3:6]) # dots of 50
	ax6 = fig.add_subplot(gs[2, 3:6]) # bar
	ax7 = fig.add_subplot(gs[3, 0:3])
	ax8 = fig.add_subplot(gs[3, 3:6])

	n = len(d_rcs)
	colors = Colors(n).colors
	for i, (sp, rcs) in enumerate(sorted(d_rcs.items())):
		_sp = convert_sp(sp)
		i50, sn50, fn50 = rcs.xn50
		fm, sm = np.median(rcs.extended_frs), np.median(rcs.extended_ns)
		om = np.median(rcs.extended_ois)
		ns = rcs.ns
		nb = len(ns) # np.sum(rcs.refcounts[:, 1]) #len(ns)
		bm = sum(rcs.refcounts[:, 1]) # sum(rcs.ns)
		kargs = dict(color=colors[i], alpha=1, label=_sp)
		ax0.plot(-2, -2, linestyle='-', marker='o',  **kargs)
#		ax0.barh(0, 0, height=0, left=0, align='center', label=_sp)

		ax1.plot(range(1, nb+1), np.cumsum(ns), drawstyle="steps-post",  **kargs)
	
		ax2.plot(range(1, nb+1), ns, drawstyle="steps-post", **kargs)
#		ax2.scatter(i50, sn50, label=_sp)
		
		ax3.scatter(len(rcs.ns), bm, **kargs)
#		ax4.hist(rcs.extended_frs, bins=40, histtype='step', **kargs)
		hist_plot(rcs.extended_frs, ax4, bins=40, **kargs)
#		ax4.scatter(rcs.frs, rcs.ns,  alpha=0.6, label=_sp)
#		ax5.scatter(fn50, sn50, alpha=0.8, label=_sp)
		ax5.scatter(fm, sm, **kargs)

#	   ax5.hist(rcs.extended_ns, bins=30, histtype='step', label=_sp)
#		ax5.scatter(np.mean(rcs.extended_frs), np.mean(rcs.extended_ns), label=_sp)
		ax6.scatter(fm, om, **kargs)
#		ax7.hist(rcs.extended_ois, bins=40, histtype='step', **kargs)
		hist_plot(rcs.extended_ois, ax7, bins=40, **kargs)

		ax8.plot(rcs.refcounts[:, 0], rcs.refcounts[:, 1], **kargs)
#		print(sn50, sm)
#	set_legends([ax1, ax2, ax3,ax4,ax5,ax6], fontsize=legend_fontsize)
#	ax0.legend(loc='upper left',fancybox=False, frameon=False, ncols=ncols)
	for ax, xlab, ylab in (
		[ax1, '# of blocks', 'Cumulative number of genes'],
		[ax2, '# of blocks', 'Number of genes'],
		[ax3, 'Total number of blocks', 'Total number of genes'],
		[ax4, 'Fractionation rate', 'Number of genes'],
		[ax5, 'Fractionation rate', 'Block size'],
		[ax6, 'Fractionation rate', 'Orthology index'],
		):
		set_labels(ax, xlab, ylab)
	
	#legend
	ncols = n // 30 + 1
	ax0.set_xlim(0,1)
	ax0.set_ylim(0,1)
	ax0.legend(loc=(-1.5, 0),fancybox=False, frameon=False, ncols=ncols)
	ax0.xaxis.set_tick_params(length=0)
	ax0.spines['right'].set_color('none')
	ax0.spines['top'].set_color('none')
	ax0.spines['left'].set_color('none')
	ax0.spines['bottom'].set_color('none')
	ax0.xaxis.set_ticks([])
	ax0.yaxis.set_ticks([])
#	ax0.set_title('Species', loc='left')

	plt.subplots_adjust(hspace=0.25, wspace=1.8)
	plt.savefig(outfig, bbox_inches='tight')

def hist_plot(data, ax, bins=10, alpha=1, **kargs):
	n,bins,patches = ax.hist(data, bins=bins, alpha=0, **kargs)
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
	return '$' + sp.replace('_', '~') + '$'
class SynRec:
	def __init__(self, rc):
		self.N = rc.N
		self.fr = rc.fr
		self.oi = rc.oi
class SynRecs:
	def __init__(self, rcs):
		self.rcs = rcs
	@lazyproperty
	def matrix(self):
		data = []
		for rc in self.rcs:
			data += [[rc.N, rc.fr, rc.oi]]
		return np.array(sorted(data, key=lambda x:-x[0]))
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
		sn50,i50 = self.get_nx0()
		fn50 = list(sorted(self.frs, reverse=1))[i50 -1]
		return i50, sn50, fn50
	@lazyproperty
	def extended_frs(self):
		frs = []
		for n, fr, oi in self.matrix:
			frs += [fr]* int(n)
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
			ois += [oi]* int(n)
		return ois

