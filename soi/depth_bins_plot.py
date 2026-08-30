import sys
import numpy as np
from math import ceil
from collections import OrderedDict
import matplotlib.pyplot as plt
from .small_tools import parse_kargs
from . import style  # global mpl settings
from .small_tools import open_file as open

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
             cmap=None, ax=None, wsize=50, figsize=5, ymax=None,
             axis='x', pos_lim=None, **kargs):
    '''Plot position vs value scatter + smoothed line.
    axis='x': position on x-axis, value on y-axis (default)
    axis='y': position on y-axis, value on x-axis (rotated)
    pos_lim: tuple (min, max) to override position axis limits
    '''
    Y1s = []
    # plot dots and smooth lines
    for x, y in zip(Xs, Ys):
        if axis == 'y':
            pos, val = list(y), list(x)
        else:
            pos, val = list(x), list(y)
        # scatter
        plt.scatter(x, y, marker=',', s=point_size,
                    c=val, cmap=cmap, alpha=alpha)
        if ymax is None:
            Y1s += val
        # smoothed median line
        if len(pos) > wsize:
            pv = [(_p, _v) for _p, _v in zip(pos, val) if _p is not None]
            pv.sort(key=lambda p: p[0])
            _pos, _val = zip(*pv)
            P_smooth, V_smooth = [], []
            for i in range(len(_pos)):
                P_smooth.append(_pos[i])
                s = int(max(0, i - wsize // 2))
                e = int(min(i + wsize // 2, len(_pos)))
                V_smooth.append(np.median(_val[s:e]))
            if axis == 'y':
                plt.plot(V_smooth, P_smooth, ls='--', lw=2, color="black")
            else:
                plt.plot(P_smooth, V_smooth, ls='--', lw=2, color="black")
    if ymax is None:
        ylim = (np.median(Y1s)+0.01) * yfold
        ymax = ylim
    # plot chrom boundary
    _lines = plt.hlines if axis == 'y' else plt.vlines
    for v in vlines:
        _lines(v, 0, ymax, color=chr_color)
    # plot chrom label (skip for axis='y' — redundant with main plot)
    if axis != 'y':
        for x, label in zip(label_x, labels):
            if vis_labels is not None and label not in vis_labels:
                continue
            y = -ymax/30
            label = label.replace('chr', '')
            plt.text(x, y, label, horizontalalignment='center',
                     verticalalignment='top', fontsize=csize)
    # custom vlines / hlines
    if vvlines is not None:
        _lines(vvlines, 0, ymax, color=arm_color,
               lw=0.5, linestyle='dashed')
    # axis labels and limits
    _val_lab = ylab if ylab is not None else 'Depth'
    if axis == 'y':
        plt.xlabel(_val_lab, fontsize=csize)
        plt.xlim(0, ymax)
        if pos_lim is not None:
            plt.ylim(pos_lim)
        else:
            plt.ylim(0, max(vlines))
        ax.yaxis.tick_right()
        ax.set_yticks([])
    else:
        plt.ylabel(_val_lab, fontsize=csize)
        plt.ylim(0, ymax)
        if pos_lim is not None:
            plt.xlim(pos_lim)
        else:
            plt.xlim(0, max(vlines))
        ax.xaxis.tick_top()
        ax.set_xticks([])
    plt.title(title, fontsize=size)
    ax.minorticks_on()
    return ax


if __name__ == '__main__':
    kargs = parse_kargs(sys.argv)
    inDepth = sys.argv[1]
    main(inDepth,  **kargs)
