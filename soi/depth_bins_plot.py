# python ~/src/depth_bins_plot3.py mapping.depth.gz.bins ../ref.fa.gap
import sys
import numpy as np
from math import ceil
from collections import OrderedDict
#from xopen import xopen as open
from .small_tools import parse_kargs

def main(inDepth=sys.argv[1], vvlines=None, order=None, window_size=50000, window_step=25000, minlength=0.005, 
		height_width_ratio=None, col=2, **kargs):
	outBinDeapth = inDepth + '.bins'
	outplot = outBinDeapth + '.pdf'
	d_bins = OrderedDict()
	d_max_pos = OrderedDict()
#	d_count = {}
	print('loading depth', file=sys.stderr)
	for line in open(inDepth):
		temp = line.rstrip().split()
		CHR = temp[0]
		POS = int(temp[1])
		DEPTH = float(temp[col])	# 2 for pos; 3 for bed
#		DEPTH = map(float, temp[3:])	# bedgraph
#		TOTAL = sum(DEPTH)
#		BIN_from = int(ceil(1.0 * (POS) / window_step))
#		BIN_to = int(ceil(1.0 * (POS+window_size) / window_step))
#		for BIN in range(BIN_from, BIN_to):
#			try: d_bins[CHR][BIN] += TOTAL
#			except KeyError: 
#				try: d_bins[CHR][BIN] = TOTAL
#				except KeyError: d_bins[CHR] = OrderedDict([(BIN, TOTAL)])
		BIN = POS
		try: d_bins[CHR][BIN] = DEPTH
		except KeyError: d_bins[CHR] = OrderedDict([(BIN, DEPTH)])
#		key = (CHR, BIN)
#		try: d_count[key] += 1
#		except KeyError: d_count[key] = 1
		d_max_pos[CHR] = max(d_max_pos.get(CHR,0), POS)

#	for k,v in d_bins.items():
#		for b, t in v.items():
#			d_bins[k][b] = 1.0 * t / d_count[(k,b)]
#	f = open(outBinDeapth, 'w')
#	for k,v in d_bins.items():
#		for b, t in v.items():
#			start = (b+1)*window_step + 1
#			print >>f, '%s\t%s\t%s' % (k, start, t)
#	f.close()
	CHRs = [line.strip().split()[0] for line in open(order)] if order is not None else list(d_bins.keys())
	CHRs = [chr for chr in CHRs if chr in d_bins]
	print(d_max_pos, file=sys.stderr)
	last_start = 0
	Xs,Ys = [], []
	labels, label_x, vlines = [], [], []
	d_offset = {}
#	for (CHR, BINs) in d_bins.items():
	for CHR in CHRs:
		BINs = d_bins[CHR]
		d_offset[CHR] = last_start
		length = d_max_pos[CHR]
		x, y = [], []
		for BIN, depth in list(BINs.items()):
#			start = (BIN+1)*window_step #window_size
			start = BIN
			start += last_start
			x += [start]
			y += [depth]
		Xs += [x]
		Ys += [y]
		last_start += length
		labels += [CHR]
		label_x += [last_start - length/2]
		vlines += [last_start]
	tot_len = sum(d_max_pos.values())
	vis_labels = set([k for k, v in list(d_max_pos.items()) if 1.0*v/tot_len >= minlength])
	if height_width_ratio is None:
		height_width_ratio = sum(d_max_pos.values()) / max(d_max_pos.values())
	if vvlines is not None:
		from minimap2synetic_plot import parse_hvlines, add_offset
		vvlines = add_offset(parse_hvlines(vvlines, min_span=200), d_offset)
		print(vvlines, len(vvlines))
	bin_plot(Xs,Ys, labels, label_x, vlines, height_width_ratio, outplot, vis_labels, vvlines=vvlines, **kargs)

def bin_plot(Xs,Ys, labels, label_x, vlines, height_width_ratio=None, 
		outplot=None, vis_labels=None, vvlines=None, ylab=None, yfold=3, 
		chr_color="black", arm_color="grey", point_size=None, title=None, 
		csize=None, size=None, alpha=None,
		cmap=None, ax=None, wsize=50, figsize=5, ymax=None, **kargs):
	import matplotlib.pyplot as plt
#	plt.figure(figsize=(5*height_width_ratio,5))
#	figsize = (figsize*height_width_ratio,figsize)
#	fig, ax = plt.subplots(1,1, figsize=figsize)
#	ny, sumy = 0, 0
	Y1s = []
#	wsize = 100	# sliding window
#	plt.scatter(Xs, Ys, marker=',', s=point_size, c = Ys, cmap = cmap, )
	for x, y in zip(Xs,Ys):
#		print(x[:10], y[:10])
#		xy = [(_x, _y) for _x, _y in zip(x,y) if isinstance(_x, int)]
#		x, y = zip(*sorted(xy))
#		plt.plot(x, y)
		plt.scatter(x, y, marker=',', s=point_size, c = y, cmap = cmap, alpha=alpha)
#		ny += len(y)
#		sumy += sum(y)
		if ymax is None:
			Y1s += y
		if len(x) > wsize:
			xy = [(_x, _y) for _x, _y in zip(x,y) if _x is not None]
			x, y = zip(*sorted(xy))

			X, Y = [], []
			for i in range(len(x)):
				X += [x[i]]
				s = int(max(0, i-wsize/2))
				e = int(min(i+wsize/2, len(x)))
#				print(s,e)
				Y += [np.median(y[s:e])] #[sum(y[s:e]) / (e-s)]
			plt.plot(X, Y, ls='--', lw=2, color="black")
#	ylim = sumy / ny * 2.5
	if ymax is None:
		ylim = (np.median(Y1s)+1) * yfold
	#	print(ylim, np.median(Y1s), yfold)
		ymax = ylim
	for v in vlines:
		plt.vlines(v, 0, ymax, color=chr_color)
	for x, label in zip(label_x, labels):
	#	if not label.startswith('chr'):
	#		continue
		if vis_labels and not label in vis_labels:
			continue
		y = -ymax/30
		label = label.replace('chr', '')
		plt.text(x, y, label, horizontalalignment='center',verticalalignment='top',fontsize=csize) #, rotation=30)

	# vlines
	if vvlines is not None:
		plt.vlines(vvlines, 0, ymax, color=arm_color, lw=0.5, linestyle='dashed')
	ylabel = ylab if ylab is not None else 'Depth'
	plt.ylabel(ylabel, fontsize=csize)
	plt.ylim(0, ymax)
	xlim = max(vlines)
	plt.xlim(0, xlim)
	ax.xaxis.tick_top()
	plt.title(title, fontsize=size)
	ax.set_xticks([])
	ax.minorticks_on()
#	plt.ticklabel_format(style='plain')
#	plt.savefig(outplot, bbox_inches='tight')

if __name__ == '__main__':
	kargs = parse_kargs(sys.argv)
	inDepth=sys.argv[1]
#	try: vlines= sys.argv[2]
#	except IndexError: vlines=None
	main(inDepth,  **kargs)
