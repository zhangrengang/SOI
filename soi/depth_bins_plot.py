import sys
import numpy as np
from math import ceil
from collections import OrderedDict
import matplotlib.pyplot as plt
# from xopen import xopen as open
from .small_tools import parse_kargs


def main(inDepth=sys.argv[1], vvlines=None, order=None,
         window_size=50000, window_step=25000, minlength=0.005,
         height_width_ratio=None, col=2, **kargs):
    outBinDeapth = inDepth + '.bins'
    outplot = outBinDeapth + '.pdf'
    d_bins = OrderedDict()
    d_max_pos = OrderedDict()
    print('loading depth', file=sys.stderr)
    for line in open(inDepth):
        temp = line.rstrip().split()
        CHR = temp[0]
        POS = int(temp[1])
        DEPTH = float(temp[col])  # 2 for pos; 3 for bed
        BIN = POS
        try:
            d_bins[CHR][BIN] = DEPTH
        except KeyError:
            d_bins[CHR] = OrderedDict([(BIN, DEPTH)])
        d_max_pos[CHR] = max(d_max_pos.get(CHR, 0), POS)

    CHRs = [line.strip().split()[0] for line in open(order)] \
        if order is not None else list(d_bins.keys())
    CHRs = [chr for chr in CHRs if chr in d_bins]
    print(d_max_pos, file=sys.stderr)
    last_start = 0
    Xs, Ys = [], []
    labels, label_x, vlines = [], [], []
    d_offset = {}
    for CHR in CHRs:
        BINs = d_bins[CHR]
        d_offset[CHR] = last_start
        length = d_max_pos[CHR]
        x, y = [], []
        for BIN, depth in list(BINs.items()):
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
    vis_labels = set([k for k, v in list(d_max_pos.items())
                      if 1.0*v/tot_len >= minlength])
    if height_width_ratio is None:
        height_width_ratio = sum(d_max_pos.values()) / max(d_max_pos.values())
    if vvlines is not None:
        from minimap2synetic_plot import parse_hvlines, add_offset
        vvlines = add_offset(parse_hvlines(vvlines, min_span=200), d_offset)
        print(vvlines, len(vvlines))
    bin_plot(Xs, Ys, labels, label_x, vlines, height_width_ratio,
             outplot, vis_labels, vvlines=vvlines, **kargs)


def bin_plot(Xs, Ys, labels, label_x, vlines, height_width_ratio=None,
             outplot=None, vis_labels=None, vvlines=None, ylab=None, yfold=2,
             chr_color="black", arm_color="grey", point_size=None, title=None,
             csize=None, size=None, alpha=None,
             cmap=None, ax=None, wsize=50, figsize=5, ymax=None, **kargs):
    Y1s = []
    # plot dots and smooth lines
    for x, y in zip(Xs, Ys):
        plt.scatter(x, y, marker=',', s=point_size,
                    c=y, cmap=cmap, alpha=alpha)
        if ymax is None:
            Y1s += y
        if len(x) > wsize:
            xy = [(_x, _y) for _x, _y in zip(x, y) if _x is not None]
            x, y = zip(*sorted(xy))

            X, Y = [], []
            for i in range(len(x)):
                X += [x[i]]
                s = int(max(0, i-wsize/2))
                e = int(min(i+wsize/2, len(x)))
                Y += [np.median(y[s:e])]
            plt.plot(X, Y, ls='--', lw=2, color="black")
    if ymax is None:
        ylim = (np.median(Y1s)+0.01) * yfold
        ymax = ylim
    # plot chrom boundary
    for v in vlines:
        plt.vlines(v, 0, ymax, color=chr_color)
    # plot chrom label
    for x, label in zip(label_x, labels):
        if vis_labels and not label in vis_labels:
            continue
        y = -ymax/30
        label = label.replace('chr', '')
        plt.text(x, y, label, horizontalalignment='center',
                 verticalalignment='top', fontsize=csize)  # , rotation=30)

    # custom vlines
    if vvlines is not None:
        plt.vlines(vvlines, 0, ymax, color=arm_color,
                   lw=0.5, linestyle='dashed')

    # ylab and ylim
    ylabel = ylab if ylab is not None else 'Depth'
    plt.ylabel(ylabel, fontsize=csize)
    plt.ylim(0, ymax)
    xlim = max(vlines)
    plt.xlim(0, xlim)
    ax.xaxis.tick_top()
    plt.title(title, fontsize=size)
    ax.set_xticks([])
    ax.minorticks_on()
    return ax


if __name__ == '__main__':
    kargs = parse_kargs(sys.argv)
    inDepth = sys.argv[1]
    main(inDepth,  **kargs)
