#!/bin/env python
#coding: utf-8
import argparse,sys, os
import re
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = 'Arial'
import logging
from .mcscan import Collinearity, Gff, XCollinearity
#from .OrthoFinder import OrthoFinder
from .ploidy_plotter import add_ploidy_opts, get_ploidy, plot_bars
from .WGDI import AK
from .RunCmdsMP import logger

__version__='0.1'
__LastModified__='20190226'
__Example__=''
def dotplot_args(parser):
	parser.add_argument('-s', metavar='FILE', type=str, required=True, nargs='+', help="syntenic block file (*.collinearity, output of MCSCANX/WGDI)[required]")
	parser.add_argument('-g', metavar='FILE', type=str, required=True, help="gene annotation gff file (*.gff, one of MCSCANX/WGDI input)[required]")
	parser.add_argument('-c', metavar='FILE', type=str, required=True, help="chromosomes config file (*.ctl, same format as MCSCANX dotplotter)[required]")
	parser.add_argument('-o', metavar='STR', type=str, default=None, help="output file prefix. [default: the same as `-c`]")
	parser.add_argument('--format', metavar='FORMAT', action='append', default=['pdf', 'png'], help="output figure format [default=%(default)s]")
	parser.add_argument('--homology', action='store_true', default=False, help="`-s` is in homology format (gene1<tab>gene2). [default=%(default)s]")
	parser.add_argument('--number-plots', action='store_true', default=False, help="number subplots with (a-d). [default=%(default)s]")
	parser.add_argument('--min-block', metavar='INT', type=int, default=None, help="min gene number in a block. [default=%(default)s]")
	parser.add_argument('--min-same-block', metavar='INT', type=int, default=25, help="min gene number in a block on the same chromosome. [default=%(default)s]")
	parser.add_argument('--min-dist', dest='tandem_dist', metavar='INT', type=int, default=None, help="remove tandem with distance shorter than this value. [default=%(default)s]")
#	parser.add_argument('--plot-dot', action='store_true', default=None, help="also plot dot without Ks. [default=%(default)s]")
#	parser.add_argument('--hide-blocks', type=str, default=None, help="blocks to hide, one block id per line. default=%(default)s")
#	parser.add_argument('--matrix', type=str, default=None, help="output chrom matrix file")
#	parser.add_argument('--plot-cluster', action='store_true', default=False, help="plot cluster tree. default=%(default)s")
#	parser.add_argument('--width', type=float, default=18, help="width of whole plot. default=%(default)s")
#	parser.add_argument('--height', type=float, default=1, help="fator of height (actual height = width*height). default=%(default)s")
#	parser.add_argument('--source', type=str, choices=['mcscanx', 'wgdi'], default=None, help="source of collinearity [default: auto]")

	group_dot = parser.add_argument_group('Dot plot', 'settings for dot plots')
	group_dot.add_argument('--cluster', action='store_true', default=False, help="cluster chromosomes. [default=%(default)s]")
	group_dot.add_argument('--diagonal', action='store_true', default=False, help="try to put blocks onto the diagonal. [default=%(default)s]")
	group_dot.add_argument('--gene-axis', action='store_true', default=False, help="use gene as axis instead of base pair. [default=%(default)s]")
	group_dot.add_argument('--xlines', metavar='FILE', type=str, default=None, help="bed/pos file to add vertical lines. [default=%(default)s]")
	group_dot.add_argument('--ylines', metavar='FILE', type=str, default=None, help="bed/pos file to add horizontal lines. [default=%(default)s]")
	group_dot.add_argument('--xbars', metavar='FILE', type=str, default=None, help="ancetor file to set colorbar for x axis. [default=%(default)s]")
	group_dot.add_argument('--ybars', metavar='FILE', type=str, default=None, help="ancetor file to set colorbar for y axis. [default=%(default)s]")
	group_dot.add_argument('--xbarlab', action='store_true', default=False, help="plot labels for x bars. [default=%(default)s]")
	group_dot.add_argument('--ybarlab', action='store_true', default=False, help="plot labels for y bars. [default=%(default)s]")

	group_dot.add_argument('--xlabel', type=str, default=None, help="x label for dot plot. [default=%(default)s]")
	group_dot.add_argument('--ylabel', type=str, default=None, help="y label for dot plot. [default=%(default)s]")
	group_dot.add_argument('--figsize', metavar='NUM', type=float, nargs='+', default=[18], help="figure size (width [height]) [default=%(default)s]")
	group_dot.add_argument('--fontsize', metavar='NUM', type=float, default=10, help="font size of chromosome labels [default=%(default)s]")
	group_dot.add_argument('--dotsize', metavar='NUM', type=float, default=1, dest='point_size', help="dot size [default=%(default)s]")


	group_orth = parser.add_argument_group('Orthology Index filter/color', 'filtering or coloring blocks by Orthology Index (prior to Ks color)')
	group_orth.add_argument('--ofdir', metavar='FOLDER/FILE', type=str, nargs='+', default=None, help="OrthoFinder output folder/ OrthoMCL output pair file. [default=%(default)s]")
	group_orth.add_argument('--of-ratio', metavar='FLOAT', type=float, default=0, help="Orthology Index cutoff [default=%(default)s]")
	group_orth.add_argument('--of-color', action='store_true', default=None, help="coloring dots by Orthology Index [default=%(default)s]")

	group_ks = parser.add_argument_group('Ks plot', 'options to histogram plot with Ks')
	group_ks.add_argument('--kaks', metavar='FILE', type=str, default=None, help="kaks output from KaKs_Calculator/WGDI. [default=%(default)s]")
	group_ks.add_argument('--ks-hist', action='store_true', default=None, help="plot histogram or not [default=%(default)s]")
	group_ks.add_argument('--max-ks', metavar='Ks', type=float, default=1, help="max Ks (x limit) [default=%(default)s]")
	group_ks.add_argument('--ks-cmap', metavar='Ks', type=float, nargs='+', default=None, help="color map for Ks. [default=%(default)s]")
	group_ks.add_argument('--ks-step', metavar='Ks', type=float, default=0.02, help="Ks step of histogram [default=%(default)s]")
	group_ks.add_argument('--use-median', action='store_true', default=False, help="use median Ks for a block. [default=%(default)s]")
	group_ks.add_argument('--method', metavar='STR', type=str, default='NG86', help='Ks calculation method [default=%(default)s]')
	group_ks.add_argument('--lower-ks', metavar='Ks', type=float, default=None, help="lower limit of median Ks. [default=%(default)s]")
	group_ks.add_argument('--upper-ks', metavar='Ks', type=float, default=None, help="upper limit of median Ks. [default=%(default)s]")
#	group_ks.add_argument('--clip-ks', action='store_true', default=None, help="clip ks > max-ks. [default=%(default)s]")
#	group_ks.add_argument('--hist-ylim', type=float, default=None, help="max y axis of Ks histgram. [default=%(default)s]")
#	group_ks.add_argument('--yn00', action='store_true', default=False, help='turn to YN00[default=%(default)s]')
#	group_ks.add_argument('--fdtv', action='store_true', default=False, help='turn to 4DTV[default=%(default)s]')

	group_ploidy = parser.add_argument_group('ploidy plot', 'options to plot relative ploidy (synteny depth)')
	group_ploidy.add_argument('--plot-ploidy', action='store_true', default=False, help="plot relative ploidy. [default=%(default)s]")
	add_ploidy_opts(group_ploidy)

	group_bin = parser.add_argument_group('plot Ks by bins', 'options to plot binned Ks')
	group_bin.add_argument('--plot-bin', action='store_true', default=False, help="plot binned Ks. [default=%(default)s]")

def reset_args(args):
	args.matrix = None
	args.hide_blocks = None
	args.plot_dot = None
	args.source = None
	args.yn00 = False
	args.fdtv = False
	args.clip_ks = None
	args.hist_ylim = None
		
	if args.o is None:
		args.o = os.path.splitext(os.path.basename(args.c))[0]
	# ploidy
	if args.window_step is None:
		args.window_step = args.window_size / 5
	if args.min_overlap is None:
		args.min_overlap = args.window_size / 2.5
	else:
		args.min_overlap = args.min_overlap*args.window_size
	if args.of_color:
		args.max_ks = min(args.max_ks, 1)
		
def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter,\
		epilog="Version: {}\nLast Modification Date: {}".format(__version__,__LastModified__),\
		version="Version: {}".format(__version__),\
		description="Example: {}".format(__Example__))

	dotplot_args(parser)
	if len(sys.argv)==1:
		args = parser.parse_args(['-h'])
	else:
		args = parser.parse_args()
	
	return args

class Args:
	def __init__(self):
		pass

def xmain(**kargs):
	args = Args()
	for k,v in kargs.items():
		setattr(args, k, v)
	return main(args)
def main(args):
	reset_args(args)
	collinearity = args.s
	gff		  = args.g
	ctl		  = args.c
	prefix	   = args.o
	kaks		 = args.kaks
#	outplot	  = prefix + '.pdf'
#	print args.format
	ks_args = {'yn00': args.yn00, 'method':args.method, 'fdtv':args.fdtv}
	if args.hide_blocks is not None:
		args.hide_blocks = set([line.strip().split()[0] for line in open(args.hide_blocks)])
	chrs1, chrs2 = parse_ctl(ctl)
	same_sp = True if chrs1 == chrs2 else False
	blocks, lines1, lines2, ortholog_graph,chrs1, chrs2, d_offset1, d_offset2 = parse_collinearity(
		collinearity, gff, chrs1, chrs2, kaks, args.homology, 
		source = args.source, 
		ofdir=args.ofdir, of_ratio=args.of_ratio, of_color=args.of_color,
		hide_blocks = args.hide_blocks, use_median=args.use_median, 
		lower_ks=args.lower_ks, upper_ks=args.upper_ks,
		cluster=args.cluster, diagonal=args.diagonal, gene_axis=args.gene_axis,
		matrix=args.matrix, min_same_block=args.min_same_block, 
		min_block=args.min_block, tandem_dist=args.tandem_dist, **ks_args)
#	logger.info('{} blocks to plot'.format(len(blocks)))	
	xpositions = [(lines1[i] + lines1[i+1]) / 2 for i in range(len(lines1)-1)]
	ypositions = [(lines2[i] + lines2[i+1]) / 2 for i in range(len(lines2)-1)]
	blocks = sorted(blocks, key=lambda x: abs(x[-1][0] - x[0][0]))	# sort by block length
	xclines = add_offset(parse_hvlines(args.xlines), d_offset1) if args.xlines else None
	yclines = add_offset(parse_hvlines(args.ylines), d_offset2) if args.ylines else None

	if kaks:
		prefix += '.ksmax' + str(args.max_ks)
		if args.use_median:
			prefix += '.median'
	if args.homology:
		prefix += '.homo'
	if args.of_ratio > 0:
		prefix += '.oimin' + str(args.of_ratio)
	if args.diagonal:
		prefix += '.diagonal'
	if kaks and args.plot_dot: # skip
		outplots = [prefix + '.dot.' + fmt for fmt in args.format]
		plot_blocks(blocks, outplots, ks = None, max_ks=None, ks_hist=None,
			xlabels=chrs1, ylabels=chrs2,
			xpositions=xpositions, ypositions=ypositions,
			xlines=lines1, ylines=lines2,
			xlim=max(lines1), ylim=max(lines2))
	ploidy_data = coord_path1, coord_path2, coord_graph1, coord_graph2 = parse_gff(gff, chrs1, chrs2)
	outplots = [prefix + '.' + fmt for fmt in args.format]
	ks = None if kaks is None and args.ofdir is None else True
	#print(lines2, yclines, d_offset2, parse_hvlines(args.ylines))
	plot_blocks(blocks, outplots, ks = ks, 
			#max_ks=args.max_ks, ks_hist=args.ks_hist, ks_cmap=args.ks_cmap, 
			#clip_ks=args.clip_ks, min_block=args.min_block, ks_step=args.ks_step,
			xlabels=chrs1, ylabels=chrs2, same_sp=same_sp,
			#hist_ylim=args.hist_ylim,
			xpositions=xpositions, ypositions=ypositions,
			xelines=lines1, yelines=lines2,
			xclines=xclines, yclines=yclines, # such as centromere
			xlim=max(lines1), ylim=max(lines2),
			xoffset=d_offset1, yoffset=d_offset2, gff=gff, 
			ploidy=args.plot_ploidy, ploidy_data = ploidy_data, ortholog_graph=ortholog_graph, **args.__dict__
			)
def is_mcscan_style(labels):
	matches = [re.compile(r'[A-Za-z]{2}\d{1,5}[A-Za-z]*$').match(lab) for lab in labels]
	return all(matches)
def _remove_prefix(labels):
	if is_mcscan_style(labels):
		return [label[2:] for label in labels]
	return labels
def plot_blocks(blocks, outplots, ks=None, max_ks=None, ks_hist=False, ks_cmap=None, clip_ks=None, min_block=None, ks_step=0.02,
			xlabels=None, ylabels=None, xpositions=None, ypositions=None, xelines=None, yelines=None, xlim=None, ylim=None,
			figsize=18, fontsize=10, point_size=0.8, xclines=None, yclines=None, plot_bin=None,
			xoffset=None, yoffset=None, xbars=None, ybars=None, gff=None, gene_axis=None, xbarlab=True, ybarlab=True, 
			hist_ylim=None, xlabel=None, ylabel=None, remove_prefix=True, number_plots=True, same_sp=False,
			ploidy=False, ploidy_data=None, ortholog_graph=None, of_color=False, homology=False, **kargs
			):
	import matplotlib
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	xcsize = ycsize = fontsize	# chromosome label
	xsize = ysize = fontsize * 2.5	# species label
	if xlabel is not None and xlabels is not None and remove_prefix:
		xlabels = _remove_prefix(xlabels)
		xcsize = xcsize*1.5 
	if ylabel is not None and ylabels is not None and remove_prefix:
		ylabels = _remove_prefix(ylabels)
		ycsize = ycsize*1.5
#	if ks_cmap:
#		cmap = create_ks_map(ks_cmap, max_ks)
#	else:
#		cmap = cm.jet
	figwidth = x = figsize[0]
	if ks is not None:
		if ks_hist is not None:
			_y = x * 1.2
			bbox_inches='standard'
		else:
			_y = x * 1.05
			bbox_inches='tight'
	else:
		_y = x * 0.98 #*1.2 if ks_hist else x
		bbox_inches='tight'
	bbox_inches='tight'
	y = _y if len(figsize) ==1 else figsize[1]
	plt.figure(figsize=(x,y))
	if ks_hist:
		ax = plt.subplot2grid((6,5),(0,0), rowspan=5, colspan=5)
	else:
		if ks is not None:
			ax = plt.subplot2grid((21,20),(0,0), rowspan=20, colspan=20)
		else:
			ax = plt.subplot2grid((5,5),(0,0), rowspan=5, colspan=5)
	# ax1
	allKs = []
	kXs, kYs, Ks = [], [], []
	for block in blocks:
		#if min_block is not None and len(block) < min_block:
		#	continue
		Xs, Ys, = [], [],
		#print(block)
		for pos1, pos2, myks in block:
			if clip_ks is not None and myks > max_ks:
				continue
			if myks is None:
				continue
			Xs += [pos1]
			Ys += [pos2]
		#	if ks <= max_ks:
		#		allKs += [ks]
		#	print(myks, max_ks)
			myks = min(myks, max_ks) #if  myks is not None else None
			Ks += [myks]
			allKs += [myks]
		kXs += Xs
		kYs += Ys
		#allKs += Ks
		if ks is None:
			plt.plot(Xs, Ys, linewidth=1.5)
		else:
			plt.plot(Xs, Ys, color="grey", ls='-', alpha=0.45, linewidth=0.55)
#		else:
			#plt.plot(Xs, Ys)
	ymin, ymax = 0, ylim
	xmin, xmax = 0, xlim

	if ks is not None:
		try:
			min_ks = min([v for v in Ks if v >0])
		except ValueError:	# ValueError: min() arg is an empty sequence
			min_ks = 0
	#	if len(Ks) > 0 and max(Ks) > max_ks:
	#		max_ks = max(Ks)
		if ks_cmap:
			cmap = create_ks_map(ks_cmap, min_ks, max_ks)
		else:
			cmap = cm.jet
	if not ks is None:
		kXs += [None, None]	# points not plot
		kYs += [None, None]
		Ks += [0, max_ks]	# unify the scale
		plt.scatter(kXs, kYs, marker=',', s=point_size, c = Ks, cmap = cmap, )
	if same_sp:
		plt.plot((xmin, xmax), (ymin, ymax), ls='--', color="grey", linewidth=0.8)

	# color bars
	xlabelpad, ylabelpad = 10, 7.5
	if xbars:
		y = ylim
		width = ylim/ 60
		ylim += width	# increase limit
		has_lab = AK(xbars).plot_dotplot(xy=y, align='edge', d_offset=xoffset, axis='x', width=width, label=xbarlab,
			gene_axis=gene_axis, gff=gff, fontsize=xcsize-1)
		if has_lab:
			xlabelpad += xcsize
	if ybars:
		x = xlim
		width = xlim/ 60
		xlim += width
		#print(yoffset)
		has_lab = AK(ybars).plot_dotplot(xy=x, align='edge', d_offset=yoffset, axis='y', width=width, label=ybarlab,
			gene_axis=gene_axis, gff=gff, fontsize=ycsize-1)
		if has_lab:
			ylabelpad += ycsize * 0.75

	# species label
	chr_color, arm_color = "dimgrey", "grey"
	if xlabel:
		ax.set_xlabel(xlabel, ha='center', fontsize=xsize, labelpad=xlabelpad)
		ax.xaxis.set_label_position('top')
	if ylabel:
		ax.set_ylabel(ylabel, rotation='vertical', ha='center', fontsize=ysize, labelpad=ylabelpad)
		ax.yaxis.set_label_position('right')

	tot_lenx, tot_leny = xlim, ylim
	# chromosome label
	for _xlabel, xposition, xline in zip(xlabels, xpositions, xelines):
		x = xline
		plt.vlines(x, ymin, ymax, color=chr_color, linewidth=1)
		x, y = xposition, -tot_leny/150.0
		plt.text(x, y, _xlabel, horizontalalignment='center',verticalalignment='top',fontsize=xcsize) #, rotation=30)
	for x in [xmin, xmax]:
		plt.vlines(x, ymin, ymax, color=chr_color, linewidth=1)
#	print(xclines, yclines)
	if xclines:
		for xline in xclines:
			plt.vlines(xline, ymin, ymax, color=arm_color, linewidth=1, ls='--')
	if yclines:
		for yline in yclines:
			plt.hlines(yline, xmin, xmax, color=arm_color, linewidth=1, ls='--')
	for _ylabel, yposition, yline in zip(ylabels, ypositions, yelines):
		y = yline
		plt.hlines(y, xmin, xmax, color=chr_color, linewidth=1)
		x, y = -tot_lenx/150.0, yposition
		plt.text(x, y, _ylabel, horizontalalignment='right',verticalalignment='center',fontsize=ycsize) #rotation=30
	for y in [ymin, ymax]:
		plt.hlines(y, xmin, xmax, color=chr_color, linewidth=1)

	plt.xlim(xmin,xlim)
	plt.ylim(ymin,ylim)
	plt.xticks([])
	plt.yticks([])
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['bottom'].set_color('none')

	lsize = 17
	if number_plots and (ks_hist or ploidy):
		label = '(a)'
		plot_label(ax, label, fontsize=lsize)

	tlabel = 'OrthoIndex' if of_color else 'Ks'
	if not ks is None and ks_hist is None and False:	# color map only
		ax = plt.subplot2grid((21,20),(20,0), colspan=5)
		plt.axis('off')
		cbar = plt.colorbar(ax=ax, orientation='horizontal', location='bottom', label=label, shrink=1, fraction=0.5)
		# cbar.ax.set_xlabel('Ks', fontsize=14)


	# ax2
	if ks_hist:
		if ploidy:
			ax = plt.subplot2grid((6,5),(5,0), colspan=3)
		else:
			ax = plt.subplot2grid((6,5),(5,0), colspan=5)
		bins = int(max_ks/ks_step)
		_xlabel = tlabel #'OrthoIndex' if of_color else 'Ks'
		#print min(allKs), max(allKs)
#		ylabel = ' of syntenic gene pairs' if homology else ' of syntenic gene pairs'
		_ylabel = ' of gene pairs'
		_histgram(ax, allKs, cmap=cmap, xlim = max_ks, ylim=hist_ylim, bins=bins, normed=False, xlabel=_xlabel, ylabel=_ylabel, fontsize=xcsize)
		label = '(b)'
		if number_plots:
			plot_label(ax, label, fontsize=lsize)
	# ax3, ax4
	if ploidy:
		coord_path1, coord_path2, coord_graph1, coord_graph2 = ploidy_data
		# bar1
		if ks_hist:
			ax = plt.subplot2grid((6,5),(5,3))
			label = '(c)'
		else:
			ax = plt.subplot2grid((6,5),(5,0), colspan=2)
			label = '(b)'
		titles = [None]
		_xlabel = 'Relative y ploidy per x'
		_ylabel = 'Number of {}-gene windows'.format(kargs['window_size'])
		plot_fold(ax, titles, coord_path1, coord_graph1, coord_graph2, ortholog_graph, xlabel=_xlabel, ylabel=_ylabel, **kargs)
		if number_plots:
			plot_label(ax, label, fontsize=lsize)
		# bar2
		if ks_hist:
			ax = plt.subplot2grid((6,5),(5,4))
			label = '(d)'
		else:
			ax = plt.subplot2grid((6,5),(5,3), colspan=2)
			label = '(c)'
		titles = [None]
		_xlabel = 'Relative x ploidy per y'
		plot_fold(ax, titles, coord_path2, coord_graph2, coord_graph1, ortholog_graph, xlabel=_xlabel, ylabel=None, **kargs)
		if number_plots:
			plot_label(ax, label, fontsize=lsize)
	plt.subplots_adjust(hspace=0.3)

	logger.info('Output figures: {}'.format(outplots))
	logging.disable()
	for outplot in outplots:
		plt.savefig(outplot, bbox_inches='tight')
	# x ~ Ks
	if plot_bin:
		outfig = os.path.splitext(outplots[0])[0] + '.bin.png'
		ymax = min(max_ks, np.percentile(Ks, 95) *1.2)
		_kargs = dict(point_size=point_size, cmap=cmap, chr_color=chr_color, arm_color=arm_color,
				ymax=ymax, figwidth=figwidth, outfig=outfig, ylab=tlabel, alpha=0.3, 
				csize=xcsize, size=xsize, wsize=kargs['window_size'])
		plot_collapse(kXs, kYs, Ks, xlabels, ylabels, xpositions, ypositions,
		    xelines, yelines, xclines, yclines, xlabel, ylabel, **_kargs)

def plot_collapse(kXs, kYs, Ks, xlabels, ylabels, xpositions, ypositions, 
		xelines, yelines, xclines, yclines, xlabel, ylabel, 
		figwidth, outfig, **kargs):
	from .depth_bins_plot import bin_plot
	height = figwidth/4
	plt.figure(figsize=(figwidth, height))
	# x
	ax = plt.subplot(2, 1, 1)
	bin_plot(Xs=[kXs], Ys=[Ks], labels=xlabels, label_x=xpositions, vlines=xelines, vvlines=xclines,
		title=xlabel, ax=ax, **kargs) 
	# y
	ax = plt.subplot(2, 1, 2)
	bin_plot(Xs=[kYs], Ys=[Ks], labels=ylabels, label_x=ypositions, vlines=yelines, vvlines=yclines,
		title=ylabel, ax=ax, **kargs)

	plt.subplots_adjust(hspace=0.6)
	plt.savefig(outfig, bbox_inches='tight')


def plot_label(ax, label, **kargs):
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	xoffset = (xmax-xmin) / 60
	yoffset = (ymax-ymin) / 110
	x = xmin - xoffset
	y = ymax + yoffset
	ax.text(x, y, label, 
		fontweight='bold',
		horizontalalignment='right',verticalalignment='bottom', **kargs)

def plot_fold(ax, titles, ref_coord_paths, ref_coord_graph, qry_coord_graph, rq_ortholog_graph, **kargs):
	d_fold = get_ploidy(ref_coord_paths, ref_coord_graph,
					qry_coord_graph, rq_ortholog_graph, **kargs)
	data = [np.array(sorted(d_fold.items()))]
	#print >>sys.stderr, kargs
	plot_bars(data, titles, ax=ax, ncol=1, nrow=1, **kargs)
def _histgram(ax, allKs, cmap=None, xlim=None, ylim=None, bins=100, normed=False, xlabel='Ks', ylabel=' of syntenic gene pairs', fontsize=10):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	import matplotlib
	if cmap is None:
		cmap = cm.jet
#	matplotlib.rcParams['xtick.minor.visible'] = True
#	matplotlib.rcParams['ytick.minor.visible'] = True
	allKs = [v for v in allKs if v >= 0 and v is not None]
	if normed:
		ylabel = 'Percent' + ylabel
	else:
		ylabel = 'Number' + ylabel
#	print allKs
	n,bins,patches = ax.hist(allKs, bins=bins, density=normed, facecolor='white', alpha=0)
	Xs, Ys = [], []
	for i in range(len(bins)-1):
		X = (bins[i] + bins[i+1])/2
		if X <= xlim:
			Xs.append((bins[i] + bins[i+1])/2)
			Ys.append(n[i])
	if ylim is None:
		ylim = 1.2*max(Ys[:-1])
	line = ax.plot(Xs, Ys, ls='--', c='grey')
#	line.set_zorder(5)
	Xs += [0, xlim] # unify the scale
	Ys += [None, None]
	point = ax.scatter(Xs, Ys, marker=',', s=14, c = Xs, cmap = cmap)
#	point.set_zorder(10)
	ax.set_xlim(0, xlim)
	ax.set_ylim(0, ylim)
	ax.set_xlabel(xlabel, fontsize=fontsize*1.1)	# Ks/OrthoIndex; fontsize
	ax.set_ylabel(ylabel, fontsize=fontsize*0.8)
	ax.minorticks_on()
	cbar = plt.colorbar(ax=ax)
	return xlim, ylim
def parse_hvlines(bedfile, min_span=10):
	lines = []
	for line in open(bedfile):
		temp = line.strip().split()
		id = temp[0]
		start = int(temp[1])
		lines += [(id, start)]
		try:
			end = int(temp[2])
			if end-start < min_span:
				continue
			lines += [(id, end)]
		except: pass
	return lines

def add_offset(positions, d_left):
	lines = []
	for id, pos in positions:
		if id not in d_left:
			logger.warn('ID: {} not in the offset'.format(id))
			continue
		lines += [d_left[id] + pos]
	return lines

def create_ks_map(ks_map, min_ks, max_ks):
	import numpy as np
	from matplotlib import cm
	from matplotlib.colors import ListedColormap, LinearSegmentedColormap
	length = 256
	#min_ks = 0 
	maps =  _norm_map(ks_map, min_ks, max_ks, length)
	print(ks_map, min_ks, max_ks, maps)
	viridis = cm.get_cmap('viridis', length)
	newcolors = viridis(np.linspace(0, 1, length))
#	tab10 = cm.get_cmap('tab10', 10).colors
	tab10 = np.array([
					  [1., 0., 1., 1.], # 深红
					  [0., 0., 1., 1.], # 蓝
					  [0., 1., 0., 1.], # 绿
#					  [0., 1., 1., 1.], # 青
#					  [1., 1., 0., 1.], # 黄
					  [1., 0., 0., 1.], # 红
#					  [1., 0., 1., 1.], # 深红
					  [1., 1., 0., 1.], # 黄
					])
	for i in range(len(maps) - 1):
		newcolors[maps[i]:maps[i+1], :] = tab10[i]
	newcmp = ListedColormap(newcolors)
	return newcmp
def _norm_map(ks_map, min_ks, max_ks, length):
	maps = [int((v-min_ks) / (max_ks-min_ks) * length)for v in ks_map if v > 0 and v < max_ks]
	if maps[0] != 0:
		maps = [0] + maps
	if maps[-1] != length:
		maps += [length]
	return maps
def parse_gff(gff, chrs1, chrs2):
	chrs = set(chrs1 + chrs2)
	coord_graph1 = nx.Graph()
	coord_graph2 = nx.Graph()
	d_gff = {}
	for line in Gff(gff):
		if not line.chrom in chrs:
			continue
		key = (line.species, line.chrom)
		try: d_gff[key] += [line]
		except KeyError: d_gff[key] = [line]

	coord_path1 = []
	coord_path2 = []
	for (sp, chrom), lines in list(d_gff.items()):
		lines = sorted(lines, key=lambda x:(x.start, -x.end))
		genes = [line.gene for line in lines]
		if chrom in set(chrs1):
			coord_path1 += [genes]
		if chrom in set(chrs2):
			coord_path2 += [genes]
		for i, line in enumerate(lines):
			if chrom in set(chrs1):
				coord_graph1.add_node(line.gene, chrom=chrom, index=i)
			if chrom in set(chrs2):
				coord_graph2.add_node(line.gene, chrom=chrom, index=i)
		for i, line in enumerate(lines[1:]):
			edge = (lines[i].gene, line.gene)
			if chrom in set(chrs1):
				coord_graph1.add_edge(*edge)
			if chrom in set(chrs2):
				coord_graph2.add_edge(*edge)
	return coord_path1, coord_path2, coord_graph1, coord_graph2

def parse_collinearity(collinearity, gff, chrs1, chrs2, kaks, homology, 
		hide_blocks=None, use_median=False, lower_ks=None, upper_ks=None,
		cluster=False, diagonal=False, gene_axis=False, source=None, 
		ofdir=None, of_ratio=0, of_color=False, tandem_dist=None, min_block=None,
		matrix=None, min_same_block=None, **ks_args):
	blocks = XCollinearity(collinearity, orthologs=ofdir, gff=gff, kaks=kaks, homology=homology, source=source, **ks_args)
#	if ofdir:
#		of = OrthoFinder(ofdir)
#		ortholog_pairs = {tuple(sorted(x)) for x in of.get_orthologs()}
#		ortholog_pairs = ortholog_pairs | {tuple(sorted(x)) for x in of.get_paralogs2()}

	chrs1s, chrs2s = set(chrs1), set(chrs2)
	d_blocks = {}
	d_blocks2 = {}
	ortholog_graph = nx.Graph()
	i,j = 0,0
	m,n = 0,0
	for rc in blocks: #.parse():
		i += 1
		m += rc.N
		if not ((rc.chr1 in chrs1s and rc.chr2 in chrs2s) or ( rc.chr1 in chrs2s and rc.chr2 in chrs1s)):
			continue
		if hide_blocks is not None and rc.Alignment in hide_blocks:
			continue
		if rc.chr1 in chrs1s and rc.chr2 in chrs2s:
			chr1, chr2 = rc.chr1, rc.chr2
			genes1, genes2 = rc.genes1, rc.genes2
			start1, start2 = rc.start1, rc.start2
		else:
			chr1, chr2 = rc.chr2, rc.chr1
			genes1, genes2 = rc.genes2, rc.genes1
			start1, start2 = rc.start2, rc.start1
		# discard some genes on the same chrom
		if not homology and min_same_block is not None and chr1==chr2 and min_same_block > rc.N:
			continue
		# discard short blocks
		if min_block is not None and rc.N < min_block:
			continue
		# discard tandem blocks
		if tandem_dist is not None and rc.is_tandem(max_dist=tandem_dist):
			continue
		# filter by Ks
		# discard median_ks < lower_ks or median_ks > upper_ks
		if lower_ks is not None and rc.median_ks < lower_ks:
			continue
		if upper_ks is not None and rc.median_ks > upper_ks:
			continue
		ks = rc.ks
		if use_median:
			ks = [rc.median_ks] * len(ks)

		if ofdir:
#			pairs = { tuple(sorted(x)) for x in rc.pairs}
#			intersect = pairs & ortholog_pairs
			ratio = rc.oi #1.0*len(intersect) / len(pairs)
			if not ratio >= of_ratio:
				continue
			if of_color:
				ks = [ratio] * len(genes1)
		ortholog_graph.add_edges_from(rc.pairs)
		try: d_blocks[(chr1, chr2)] += [(genes1, genes2, ks)]
		except KeyError: d_blocks[(chr1, chr2)] = [(genes1, genes2, ks)]
		if chr1 in chrs2s and chr2 in chrs1s:	# same species
			try: d_blocks[(chr2, chr1)] += [(genes2, genes1, ks)]
			except KeyError: d_blocks[(chr2, chr1)] = [(genes2, genes1, ks)]
		if diagonal:
			value = [rc.score, start1, start2, rc.median_ks]
			try: d_blocks2[(chr1, chr2)] += [value]
			except KeyError: d_blocks2[(chr1, chr2)] = [value]
		j += 1
		n += rc.N
	logger.info('retained: {}/{} blocks, {}/{} genes'.format(j, i, n,m))
	
	if diagonal:
		chrs1, chrs2 = diagonal_chroms(d_blocks2, chrs1, chrs2)
	if cluster:
		chrs1, chrs2 = cluster_chroms(d_blocks, chrs1, chrs2)
	
	d_length = rc.chr_length
	if gene_axis:
		d_length = rc.chr_ngenes
	d_offset1, lines1 = _offset(chrs1, d_length)
	d_offset2, lines2 = _offset(chrs2, d_length)

	xblocks = []
	ksx = []
	for chr_pair, tblocks in list(d_blocks.items()):
		chr1, chr2 = chr_pair
		for (genes1, genes2, ksS) in tblocks:
			points = []
			for gene1, gene2, ks in zip(genes1, genes2, ksS):
#				if ks == 0:
#					print gene1, gene2
				if isinstance(gene1, str) or isinstance(gene2, str):
					continue
				if gene_axis:
#					print [gene1], gene1.index, d_offset1[chr1]
					pos1 = gene1.index + 1 + d_offset1[chr1]
					pos2 = gene2.index + 1 + d_offset2[chr2]
				else:
					pos1 = (gene1.start + gene1.end) / 2 + d_offset1[chr1]
					pos2 = (gene2.start + gene2.end) / 2 + d_offset2[chr2]
				points += [(pos1, pos2, ks)]
				ksx += [ks]
			if not points:
				continue
			xblocks += [points]
#	print(xblocks)
	if len(xblocks) == 0:
		logger.warn('No genes are retained or can be found in `Gff`. Check your files.')
	ksx = set(ksx)
	if len(ksx) == 1:
		logger.warn('All Ks have the same value: {}. Ks color map will be disabled'.format(ksx))
	if matrix is not None:
		fout = open(matrix, 'w')
		d_matrix = to_matrix(d_blocks)
		line = ['#chrom', 'gene'] + chrs2
		print('\t'.join(line), file=fout)
		for chr1 in chrs1:
			for gene1 in blocks.d_chrom[chr1]:
				genes2 = [d_matrix.get((chr1, chr2, gene1), '')for chr2 in chrs2]
				line = [chr1, gene1.id] + genes2
				print('\t'.join(line), file=fout)
		fout.close()
	return xblocks, lines1, lines2, ortholog_graph, chrs1, chrs2, d_offset1, d_offset2
	
def diagonal_chroms(d_blocks, chrs1, chrs2, **kargs):
	d_distance = {}
	
	for (chr1, chr2), values in list(d_blocks.items()): #[rc.score, rc.start1, rc.start2, rc.median_ks]
		score = sum([value[0] for value in values])
		start1 = min([value[1] for value in values])
		start2 = min([value[2] for value in values])
		ks = np.median([value[2] for value in values])
		d_distance[(chr1, chr2)] = (score, start1, ks)
		d_distance[(chr2, chr1)] = (score, start2, ks)
	if len(chrs1) > len(chrs2):	# 1 diagonal
		#chrs2 = best_match(d_distance, chrs1, chrs2)
		chrs1 = best_match(d_distance, chrs2, chrs1)
	else:	# 2 diagonal
		#chrs1 = best_match(d_distance, chrs2, chrs1)
		chrs2 = best_match(d_distance, chrs1, chrs2)
	return chrs1, chrs2
	
def best_match(d_distance, chrs1, chrs2): 	# ref, qry
	retain = []
	outside = []
	# best match = max score
	for chr2 in chrs2:
		blocks = [(chr1, d_distance.get((chr1, chr2), (0,0,0)) ) for chr1 in chrs1]
		best = max(blocks, key=lambda x:x[1][0])
		if best[1][0] == 0:
			outside += [chr2]
		else:
			retain += [(chr2,) + best]
	# sort by chr1 and coord
	retain = sorted(retain, key=lambda x:(chrs1.index(x[1]), x[2][1]))
	chrs2 = [val[0] for val in retain] + outside
	return chrs2
	
def cluster_chroms(d_blocks, chrs1, chrs2, ax=None, **kargs):
	d_distance = {}
	for (chr1, chr2), values in list(d_blocks.items()):
		ks = [value[2] for value in values]
		distance = calculate_distance(ks, **kargs)
		d_distance[(chr1, chr2)] = d_distance[(chr2, chr1)] = distance
	distances1 = get_distance(d_distance, chrs1, chrs2)
	distances2 = get_distance(d_distance, chrs2, chrs1)
	P1 = hierarchy(distances1)
	P2 = hierarchy(distances2)
	ordered_chrs1 = [chrs1[idx] for idx in P1['leaves']]
	ordered_chrs2 = [chrs2[idx] for idx in P2['leaves']]
	return ordered_chrs1, ordered_chrs2
def get_distance(d_distance, chrs1, chrs2):
	distances = []
	for chr1 in chrs1:
		distances += [[d_distance.get((chr1, chr2),1) for chr2 in chrs2]]
	return distances
def hierarchy(X):
	import scipy.cluster.hierarchy as sch
	X = np.array(X)
	d = sch.distance.pdist(X)
	Z = sch.linkage(d,method='average')
	P = sch.dendrogram(Z)
	return P

def calculate_distance(xks, min_ks=0, max_ks=3, **kargs):
	distance = 1.0
	n = 0
#	print(xks[:10])
	for ks in xks:
#		print(min_ks, ks, max_ks)
		try: ks = np.median(ks)
		except TypeError: continue
		if not min_ks < ks < max_ks:
			continue
		distance = distance * ks
		n += 1
	if n > 0:
		distance = distance**(1.0/n) / (n**0.5)
	else:
		distance = max_ks
	return distance
	
def to_matrix(d_blocks, ):
	d_genes = {}
	for (chr1, chr2), blocks in list(d_blocks.items()):
		blocks = sorted(blocks, key=lambda x:len(x[0]))
		for (genes1, genes2, _) in blocks:
			for gene1, gene2 in zip(genes1, genes2):
				d_genes[(chr1, chr2, gene1)] = gene2.id
	return d_genes

def _offset(chrs, d_length):
	d = {}
	last = 0
	lines = [last]
	for chr in chrs:
		offset = last
		d[chr] = offset
		last = offset + d_length[chr]
		lines += [last]
	return d, lines
def parse_ctl(ctl):
	i = 0
	chrs2 = []
	for line in open(ctl):
		line = line.split('//')[0]
		if not line:
			continue
		i += 1
		if i == 3: # x
			chrs1 = line.strip().strip(',').split(',')
			chrs1 = list(map(strip_blank, chrs1))
		if i >= 4: # y
			_chrs2 = line.strip().strip(',').split(',')
			chrs2 += list(map(strip_blank, _chrs2))
	if i < 4:
		raise ValueError('*.ctl file has without >= 4 lines')
	#chrs1, chrs2 = map(lambda x:x.strip(), [chrs1, chrs2])
#	return chrs1, chrs2
	return chrs2, chrs1
def strip_blank(x):
	return x.strip()
if __name__ =='__main__':
	main(makeArgparse())
	
