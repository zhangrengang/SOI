import sys
import os
import re
import itertools
import numpy as np
from collections import Counter, OrderedDict
from Bio import SeqIO
import matplotlib.pyplot as plt
from .mcscan import Gff, XGff, ColinearGroups, KaKsParser, seq2dict, Collinearity
from .OrthoFinder import catAln, format_id_for_iqtree
from lazy_property import LazyWritableProperty as lazyproperty
from .small_tools import mkdirs, flatten, test_s, test_f, parse_kargs
from .RunCmdsMP import run_cmd, run_job, logger
from .style import format_species_label, format_species_title

sg_colors = ['#f9c00c', '#00b9f1', '#7200da', '#f9320c', '#00b8a9']


class TPM:
    '''Parser for TMP table'''

    def __init__(self, tpmfile):
        self.tpmfile = tpmfile

    def __iter__(self):
        return self._parse()

    def _parse(self):
        for i, line in enumerate(open(self.tpmfile)):
            if i == 0 or line.startswith('gene'):
                continue
            yield TPMLine(line)

    def to_dict(self):
        d = {}
        for line in self:
            d[line.gene] = line.mean
        return d


class TPMLine:
    def __init__(self, line):
        temp = line.strip().split()
        self.gene = temp[0]
        self.exps = list(map(float, temp[1:]))

    @property
    def mean(self):
        return np.mean(self.exps)


class AK:
    """Parser for WGDI ancestor.txt+"""

    def __init__(self, akfile):
        self.akfile = akfile

    def __iter__(self):
        return self._parse()

    def _parse(self):
        segments = []
        for line in self._iter_lines():
            segment = Segment(line)
            if segments and segment.chrom != segments[-1].chrom:
                yield Chromosome(segments)
                segments = []
            segments += [segment]
        yield Chromosome(segments)

    def _parse_line(self):
        for line in self._iter_lines():
            yield Segment(line)

    def _iter_lines(self):
        for line in open(self.akfile):
            if line.startswith('#'):
                continue
            yield line

    def map_to_raw_coord(self, gff, chrmap=None):
        d_chr = Collinearity(chrmap=chrmap).map_chr()
        d_index = XGff(gff).get_index()
        for line in self._parse_line():
            gene1 = d_index[(line.chrom, line.start)]
            gene2 = d_index[(line.chrom, line.end)]
            start, end = gene1.start, gene2.end
            chrom = d_chr.get(line.chrom, line.chrom)
            line.chrom, line.start, line.end = chrom, start, end
            yield line

    def lazy_lines(self, gene_axis=True, **kargs):
        is_order = self.is_order()
        if not gene_axis and is_order:
            return self.map_to_raw_coord(**kargs)
        return self._parse_line()

    def is_order(self, ):
        ends = [seg.end for seg in self._parse_line()]
        if max(ends) > 2e6:
            return False
        return True

    def plot(self, sort_by_sg=True, outfig=None, bar_width=0.5,
             hatchs=['//', r'\\', '++', 'xx', '--'],
             plot_labels=True, min_segment=50):
        chromosomes = list(self)
        subgenomes = set([])
        for chromosome in chromosomes:
            for segment in chromosome:
                subgenomes.add(segment.subgenome)
        subgenomes = sorted(subgenomes)
        if len(subgenomes) == 1 or not hatchs:
            hatchs = [None] * len(subgenomes)
        if len(subgenomes) == 1:
            sort_by_sg = False
        d_hatchs = dict(list(zip(subgenomes, hatchs)))
        chromosomes = sorted(chromosomes, key=lambda x: int(x.label))
        if sort_by_sg:
            chromosomes = sorted(
                chromosomes, key=lambda x: x.longest_subgenome)
        labels = [chrom.label for chrom in chromosomes]
        if not outfig:
            outfig = self.akfile + '.pdf'

        import matplotlib.pyplot as plt
        import matplotlib as mpl

        Nchrom = len(chromosomes)
        width = Nchrom*1.0 / 3
        plt.figure(figsize=(width, 5))
        ax = plt.subplot(111)
        for i, chrom in enumerate(chromosomes):
            label = chrom.chrom
            yoffset = 0
            for seg in chrom:
                color, length, hatch = seg.color, len(
                    seg), d_hatchs[seg.subgenome]
                if length < min_segment:
                    color = 'white'
                plt.bar(i, length, width=bar_width, color=color,
                        bottom=yoffset, hatch=hatch, align='center')
                yoffset += length
        ax.tick_params('x', which='major', length=0, labelsize=15)
        ax.xaxis.set_tick_params(length=0)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        xlabel = '$x$ = {}'.format(Nchrom)
        plt.ylabel('')
        plt.savefig(outfig, bbox_inches='tight', transparent=True)

    def plot_dotplot(self, ax=None, d_offset={}, gff=None, align='center', xy=1,  axis='x', width=0.5, label=True, gene_axis=False, fontsize=10, color_fn=None):
        bar = plt.bar if axis == 'y' else plt.barh
        has_lab = False
        texts = []
        for sgement in self.lazy_lines(gene_axis, gff=gff):
            if sgement.chrom not in d_offset:
                continue
            offset = d_offset[sgement.chrom] + sgement.start
            lab = sgement.label if label else None
            color = color_fn(sgement) if color_fn else sgement.color
            bar(xy, len(sgement), width, offset,
                color=color, align=align, )
            if not lab:
                continue
            has_lab = True
            if axis == 'x':
                x = offset + len(sgement)/2
                y = xy + width*1.05
                text = plt.text(x, y, lab, horizontalalignment='center',
                                verticalalignment='bottom', fontsize=fontsize)
            elif axis == 'y':
                x = xy + width*1.07
                y = offset + len(sgement)/2
                text = plt.text(x, y, lab, horizontalalignment='left',
                                verticalalignment='center', fontsize=fontsize, rotation=90)
            texts += [text]
        return has_lab

    def to_jcvi(self, gene_axis=False, gff=None, outpre='jcvi'):
        f1 = open(outpre+'.bed', 'w')
        f2 = open(outpre+'.idmap', 'w')
        for seg in self.lazy_lines(gene_axis, gff=gff):
            line = [seg.chrom, seg.start, seg.end, seg.id]
            print('\t'.join(map(str, line)), file=f1)
            line = [seg.id, seg.id, seg.color]
            print('\t'.join(map(str, line)), file=f2)
        f1.close()
        f2.close()


def _adjust_text(*args):
    from adjustText import adjust_text
    adjust_text(*args)


class Chromosome:
    def __init__(self, segments):
        self.segments = segments

    def __iter__(self):
        return iter(self.segments)

    def __len__(self):
        return sum(map(len, self.segments))

    @property
    def chrom(self):
        return self.segments[0].chrom

    @property
    def label(self):
        return re.compile(r'(\d+)').search(self.chrom).groups()[0]

    @property
    def longest_subgenome(self):
        super_segments = []
        for subgenome, segments in itertools.groupby(
                sorted(self.segments, key=lambda x: x.subgenome),
                key=lambda x: x.subgenome):
            segments = Chromosome(list(segments))
            segments.subgenome = subgenome
            super_segments += [segments]
        return max(super_segments, key=lambda x: len(x)).subgenome


class Segment:
    def __init__(self, line):
        temp = line.strip().split()
        self.chrom, self.start, self.end, self.color, self.subgenome = \
            temp[:5]
        self.start, self.end = int(self.start), int(self.end)
        self.subgenome = int(self.subgenome)
        try:  #plus format
            self.label = temp[5]
        except IndexError:
            self.label = None

    def __len__(self):
        return self.end - self.start + 1

    @lazyproperty
    def id(self):
        return '{}:{}-{}|{}'.format(self.chrom, self.start, self.end, self.subgenome)

    def __repr__(self):
        return self.id

    def write(self, fout):
        print('\t'.join(map(str, (self.chrom, self.start, self.end, self.color,
                        self.subgenome))), file=fout)


def stats_besthits(alignment, pol_indice, dip_indice, labels, blasts):
    from Blast import MultiBlastOut
    logger.info('loading {}'.format(blasts))
    d_besthits = MultiBlastOut(blasts).filter_besthit(pair=True)
    logger.info('loading {}'.format(alignment))
    polAlns = Alignment(alignment, indice=pol_indice)
    dipAlns = Alignment(alignment, indice=dip_indice)
    d_count = {}
    x, y = 0, 0
    for polAln, dipAln in zip(polAlns, dipAlns):
        if dipAln.missing_rate > 0:
            continue
        i = 0
        for polG, polSP in zip(polAln.genes, polAln.species):
            i += 1
            if polG is None:
                continue
            y += 1
            bls = []
            j = 0
            for dipG, dipSP in zip(dipAln.genes, dipAln.species):
                j += 1
                key = (polG, dipG)
                if not key in d_besthits:
                    continue
                bl = d_besthits[key]
                bl.sp, bl.idx = dipSP, j
                bls += [bl]
            if len(bls) != j:
                x += 1
                if x < 10:
                    print(key, file=sys.stderr)
                elif x == 10:
                    print('...', file=sys.stderr)
                continue
            if len(bls) == 0:
                print(key, file=sys.stderr)
                continue
            best = max(bls, key=lambda x: x.bitscore)
            bestSp = best.sp
            key1, key2 = (i, polSP), (best.idx, bestSp)
            try:
                d_count[key1][key2] += 1
            except KeyError:
                try:
                    d_count[key1][key2] = 1
                except KeyError:
                    d_count[key1] = {key2: 1}
    print(x, 'errors;', y, 'valid', file=sys.stderr)
    line = ['sg_i', 'sg_lab', 'polyploid',
            'prog_i', 'diploid', 'count', 'percent']
    print('\t'.join(map(str, line)))
    for (i, polSP), d_count2 in sorted(d_count.items()):
        s_count = 0
        for (j, bestSp), count in sorted(d_count2.items()):
            s_count += count
        for (j, bestSp), count in sorted(d_count2.items()):
            line = i, labels[i-1], polSP, j, bestSp, count, 1e2*count/s_count
            print('\t'.join(map(str, line)))


def set_color_for_jcvi(alignment, bed,  gff, simples):
    d_genes = Gff(gff).get_indexed_genes()
    alns = Alignment(alignment)
    for line in open(bed):
        temp = line.strip().split()
        sp, g1, g2, color = temp[:4]
        lines = alns.fetch(g1, g2)
        array = np.array([line.genes for line in lines])
        blocks = []
        for col in range(array.shape[1]):
            genes = array[:, col]
            block = group_genes(genes, d_genes)
            block = (block[0], block[-1])
            print(block)
            blocks += [block]
        for simple in simples:
            outfile = simple + '.2'
            overlap_and_color_simple(simple, outfile, blocks, d_genes, color)


def overlap_and_color_simple(simple, outfile, blocks, d_genes, color):
    f = open(outfile, 'w')
    for line in open(simple):
        line = line.split('*')[-1]
        temp = line.strip().split()
        b1g1, b1g2, b2g1, b2g2 = temp[:4]
        if overlap_blocks(b1g1, b1g2, blocks, d_genes) and \
                overlap_blocks(b2g1, b2g2, blocks, d_genes):
            line = color + '*' + line
        f.write(line)
    f.close()


def overlap_blocks(g1, g2, blocks, d_genes, min_cov=0.1):
    g1, g2 = d_genes[g1], d_genes[g2]
    if g1.index > g2.index:
        g1, g2 = g2, g1
    for xg1, xg2 in blocks:
        xg1, xg2 = d_genes[xg1], d_genes[xg2]
        if len({g1.chrom, g2.chrom, xg1.chrom, xg2.chrom}) > 1:
            continue
        if xg1.index > xg2.index:
            xg1, xg2 = xg2, xg1
        cov = min(g2.index, xg2.index) - max(g1.index, xg1.index)
        if not cov > 0:
            continue
        if 1.0 * cov / abs(g2.index-g1.index) < min_cov and \
                1.0 * cov / abs(xg2.index-xg1.index) < min_cov:
            continue
        return True
    return False


def group_genes(genes, d_genes, max_dist=50):
    genes = [g for g in genes if g]
    genes = sorted(genes, key=lambda x: (d_genes[x].chrom, d_genes[x].index))
    gene = genes[0]
    last_chrom, last_index = d_genes[gene].chrom, d_genes[gene].index
    bins = [[gene]]
    for gene in genes[1:]:
        chrom, index = d_genes[gene].chrom, d_genes[gene].index
        if chrom == last_chrom and index-last_index < max_dist:
            bins[-1] += [gene]
        else:
            bins += [[gene]]
        last_chrom, last_index = chrom, index
    return max(bins, key=lambda x: len(x))



def plot_depth_violin(depthfile, outfigs=None, retained_only=True, drop_zero=True,
                      normalize=True, norm_by='mean', ylim=(0, 2.5), colors=None,
                      alpha=0.5, ref=None, sample=None, figsize=None, dpi=300, **kargs):
    '''每个亚基因组拷贝的深度分布：小提琴图 + 中值线（compmap 深度表专用）。

    depthfile: soi-wgdi compmap 产出的 <prefix>.depth.tsv，列 gene/copy/retained/depth。
        新版（流式）表带表头且只含完全保留行；旧版（非流式）表无表头、含全部参考基因
        （非保留行 copy=NA、retained=0），两种都吃。
    默认只取 retained==1（完全保留）的位点，并先标准化：深度 / 所选位点深度的总平均。
    于是各拷贝可直接横向比较：中值齐平在 1 附近 = 深度均匀（单向判据，不均匀 -> 不共享
    WGD，反之不成立）。
    outfigs 给 .pdf 路径，同时出同名 .png；返回 [(copy, 中值), ...]。
    '''
    d_copy = {}
    nzero = 0
    for line in open(depthfile):
        arr = line.rstrip('\n').split('\t')
        if len(arr) < 4:
            continue
        gene, copy, retained, depth = arr[:4]
        if gene == 'gene' and copy == 'copy':       # 表头
            continue
        if retained_only and retained != '1':
            continue
        if copy in ('', 'NA'):
            continue
        try:
            d = float(depth)
        except ValueError:
            continue
        if drop_zero and d == 0:      # 深度 0 = 该位点在这份 reads 里没覆盖，画分布要剔掉
            nzero += 1
            continue
        d_copy.setdefault(copy, []).append(d)
    if not d_copy:
        raise ValueError('no retained depth records in {}'.format(depthfile))
    if nzero:
        logger.info('depth violin: dropped {} zero-depth loci'.format(nzero))
    copies = list(d_copy)        # 保持表内出现顺序（= rename.list 的拷贝顺序）
    d_arr = {c: np.array(d_copy[c], dtype=float) for c in copies}
    nloci = min(len(d_arr[c]) for c in copies)
    # 深度分布有长尾（本例最大 131.8x）：y 轴上限定死，下界按分位取（数据驱动）
    pooled = np.concatenate([d_arr[c] for c in copies])
    if normalize:
        if norm_by == 'median':          # 长尾（靶向捕获效率差异大）时用中位数当分母更稳
            mu = float(np.median(pooled))
        elif norm_by == 'mean':
            mu = float(pooled.mean())
        else:
            raise ValueError("norm_by must be 'mean' or 'median', got {}".format(norm_by))
    else:
        mu = 1.0
    lo = 0.0 if (not ylim or ylim[0] is None) else float(ylim[0])
    hi = float(pooled.max()) / mu if (not ylim or ylim[1] is None) else float(ylim[1])
    if figsize is None:
        # 宽度随亚基因组数量调节（每个拷贝 0.75 英寸 + 两侧余量），整体比原来小
        figsize = (0.75 * len(copies) + 1.3, 3.6)
    if colors is None:
        colors = sg_colors
    # 自己算 KDE 并用 fill_betweenx 画：matplotlib 的 violinplot 把密度网格铺在
    # 数据全量程上（本例 0-131x），可见窗口里只剩两三个采样点 -> 平滑曲线变折线；
    # 这里把网格固定铺在 [lo, hi]（256 点），数据仍用全量（不在边界堆积）
    from scipy.stats import gaussian_kde
    ygrid = np.linspace(lo, hi, 256)
    fig, ax = plt.subplots(figsize=figsize)
    for i, c in enumerate(copies):
        v = d_arr[c] / mu
        dens = gaussian_kde(v)(ygrid)
        w = dens / dens.max() * 0.4        # 每把琴按自身峰值定宽
        ax.fill_betweenx(ygrid, i + 1 - w, i + 1 + w,
                         facecolor=colors[i % len(colors)], alpha=alpha,
                         edgecolor='#333333', linewidth=0.6, zorder=2)
        m = float(np.median(v))
        ax.plot([i + 1], [m], marker='o', color='black', markersize=4.2,
                linestyle='none', zorder=4)
        logger.info('copy {}: n={} mean {:.2f}x median {:.3f}x'.format(
            c, len(d_arr[c]), d_arr[c].mean() / mu, m))
    ax.set_ylim(lo, hi)
    logger.info('depth violin: {} x {} loci/copy, norm_by={} (divisor {:.3f}), ylim {:.3f}-{:.3f}'.format(
        len(copies), nloci, norm_by, mu, lo, hi))
    ax.axhline(1.0, color='#888888', ls='--', lw=0.8, zorder=1)
    ax.set_xticks(range(1, len(copies) + 1))
    ax.set_xticklabels(copies)
    ax.set_xlabel('Subgenome')
    ax.set_ylabel('Normalized depth' if normalize else 'Depth')
    sub = 'reads competitively mapped onto the reference subgenomes'
    if normalize:
        sub += '; depth normalized to the {}'.format(norm_by)
    if ref and sample:
        # 左 = 测序物种（reads），右 = 参考物种（其亚基因组拷贝）；用 / 分隔，整串可斜体
        t, tstyle = format_species_title('{} / {}'.format(sample, ref), sep=' / ')
        ax.set_title(t, fontstyle=tstyle, pad=20)
#        ax.text(0.5, 1.01, sub, transform=ax.transAxes,
#                ha='center', va='bottom', fontsize=9)
    elif ref:
        t, tstyle = format_species_label(ref)
        ax.set_title(t, fontstyle=tstyle, pad=20)
        ax.text(0.5, 1.01, sub, transform=ax.transAxes,
                ha='center', va='bottom', fontsize=9)
    elif sample:
        ax.set_title(sample, pad=20)
    if outfigs:
        fig.savefig(outfigs, bbox_inches='tight')
        fig.savefig(str(outfigs)[:-4] + '.png', dpi=dpi, bbox_inches='tight')
        logger.info('Written {} (+ .png)'.format(outfigs))
    plt.close(fig)
    return [(c, float(np.median(d_arr[c] / mu))) for c in copies]

class Alignment:
    '''Parser for WGDI -a output (alignment file)'''

    def __init__(self, alignment, indice=None, idmap=None, colnames=None, **kargs):
        self.alignment = alignment
        self.indice = self.parse_indice(indice)  # cut columns
        if idmap is None and colnames is not None:
            idmap = colnames
        try:
            self.idmap = self.parse_idmap(idmap)  # idx - name * unique ID
        except:
            self.idmap = None
        self.colnames = self.parse_colnames(colnames)  # names    * repeat allowed
        if self.indice is not None and colnames is not None:
            # cut colnames along with indice, keep consistent with line.genes
            self.colnames = [self.colnames[i] for i in sorted(self.indice)]

    def __iter__(self):
        return self._parse()

    def indexed(self, start=1):
        '''Iterate lines together with their row number (1-based). This is the
        numbering behind the 'XG' ids shared by compmap and genetrees, so ids
        of the same row agree across subcommands.'''
        for i, line in enumerate(self._parse(), start):
            yield i, line

    def _parse(self):
        for line in open(self.alignment):
            yield AlignmentLine(line, self.indice)

    def parse_colnames(self, colnames):
        if colnames is None:
            return self.species
        names = []
        for line in open(colnames):
            tmp = line.strip().split()
            try:
                name = tmp[1]
            except IndexError:
                name = tmp[0]
            names += [name]
        return names

    def to_array(self):
        return np.array([line.genes for line in self])

    def stat_sps(self, prefix='stat_sps'):  # by species
        f1 = open(prefix+'.copy', 'w')
        f2 = open(prefix+'.ngenes', 'w')
        line = ['copy', 'freq', 'column']
        print('\t'.join(line), file=f1)
        line = ['freq', 'column']  # total gene number, sp
        print('\t'.join(line), file=f2)

        for sp in sorted(set(self.colnames), key=lambda x: self.colnames.index(x)):
            idx = [i for i, nm in enumerate(self.colnames) if nm == sp]
            self.indice = idx
            n, copy = 0, []
            for line in self:
                genes = [g for g in line.genes if g is not None]
                ng = len(genes)
                n += ng
                copy += [ng]
            for c, freq in list(Counter(copy).items()):
                line = [c, freq, sp]
                line = list(map(str, line))
                print('\t'.join(line), file=f1)
            line = [n, sp]
            line = list(map(str, line))
            print('\t'.join(line), file=f2)

    def parse_by_species(self):
        for sp in sorted(set(self.colnames), key=lambda x: self.colnames.index(x)):
            idx = [i for i, nm in enumerate(self.colnames) if nm == sp]
            self.indice = idx
            yield sp, self, idx

    def get_high_retain(self, min_retain=1, fout=sys.stdout, **kargs):
        line = ['gene', 'species', 'group', 'xg']
        print('\t'.join(line), file=fout)
        for sp, lines, _ in self.parse_by_species():
            i = 0
            for line in lines:
                i += 1
                if 1 - line.missing_rate < min_retain:
                    continue
                for gene, xsp in line.nonmissing:
                    gname = gene.split('|', 1)[1]
                    line = [gname, xsp, sp, i]
                    line = list(map(str, line))
                    print('\t'.join(line), file=fout)

    def get_high_loss(self, max_copy=1, fout=sys.stdout, **kargs):
        line = ['gene', 'species', 'group', 'xg']
        print('\t'.join(line), file=fout)
        for sp, lines, _ in self.parse_by_species():
            i = 0
            for line in lines:
                i += 1
                genes = [g for g in line.genes if g is not None]
                if len(genes) > max_copy:
                    continue
                for gene, xsp in line.nonmissing:
                    gname = gene.split('|', 1)[1]
                    line = [gname, xsp, sp, i]
                    line = list(map(str, line))
                    print('\t'.join(line), file=fout)

    def to_bimat(self, fout=sys.stdout, start=1):
        array = self.to_array()
        for col in range(start, array.shape[1]):
            genes = array[:, col]
            id = self.colnames[col]
            seq = []
            for gene in genes:
                bs = 0 if gene is None else 1
                seq += [bs]
            seq = ''.join(map(str, seq))
            print('>{}\n{}'.format(id, seq), file=fout)

    def stat_sgs(self, prefix='stat_sgs'):  # by subgenome
        array = self.to_array()
        f1 = open(prefix+'.loss', 'w')
        f2 = open(prefix+'.retain', 'w')
        f3 = open(prefix+'.retained.all', 'w')

        line = ['loss', 'freq', 'column']
        print('\t'.join(line), file=f1)
        line = ['retain', 'freq', 'column']
        print('\t'.join(line), file=f2)
        line = ['freq', 'column']
        print('\t'.join(line), file=f3)

        for col in range(array.shape[1]):
            genes = array[:, col]
            loss, retain = [0], [0]
            xfreq = 0
            for gene in genes:
                if gene is None:
                    loss[-1] += 1
                    retain += [0]
                else:
                    xfreq += 1
                    loss += [0]
                    retain[-1] += 1
            for c, freq in list(Counter(loss).items()):
                line = [c, freq, self.colnames[col]]
                line = list(map(str, line))
                print('\t'.join(line), file=f1)
            for c, freq in list(Counter(retain).items()):
                line = [c, freq, self.colnames[col]]
                line = list(map(str, line))
                print('\t'.join(line), file=f2)
            line = [xfreq, self.colnames[col]]
            line = list(map(str, line))
            print('\t'.join(line), file=f3)
        for f in (f1, f2, f3):
            f.close()

    def parse_idmap(self, idmap):
        if idmap is None:
            return idmap
        ids = [line.strip().split()[1] for line in open(idmap)]
        if self.indice is not None:
            ids = [ids[i] for i in sorted(self.indice)]
        return ids

    def parse_indice(self, indice):
        if indice is None:
            return indice
        i = []
        for x in indice.strip(',').split(','):
            x = x.strip()
            if '-' in x:
                a, b = x.split('-')
                i += list(range(int(a)-1, int(b)))
            else:
                i += [int(x)-1]
        return i

    @property
    def species(self):
        sps = []
        for line in self:
            if not sps:
                sps = line.species
                continue
            idx = [i for i, sp in enumerate(sps) if sp is None]
            if not idx:
                break
            for i in idx:
                sps[i] = line.species[i]
        # number
        d_count = Counter(sps)
        d_counter = {sp: 0 for sp in d_count}
        for i, sp in enumerate(sps):
            if d_count[sp] == 1:
                continue
            d_counter[sp] += 1
            sp = '{}_{}'.format(sp, d_counter[sp])
            sps[i] = sp
        return sps

    def sg_exp(self, expfile, norm=False, max_missing=0, fout=sys.stdout, min_tpm=1, **kargs):
        d_exp = TPM(expfile).to_dict()
        logger.info('{} TPM loaded'.format(len(d_exp)))
        for sp, lines, idx in self.parse_by_species():
            for line in lines:
                if line.missing_rate > max_missing:
                    continue
                genes = [g if g is None else g.split(
                    '|', 1)[-1] for g in line.genes]
                exps = [d_exp.get(g, 0) for g in genes]
                _sum = sum(exps)
                if _sum < min_tpm:
                    continue
                if norm:  # normalize
                    exps = [1.0*e/_sum for e in exps]
                for i, exp in enumerate(exps):
                    if not norm and exp == 0:
                        continue
                    try:
                        sg = self.idmap[idx[i]]
                    except:
                        sg = 'sg{}'.format(i)
                    line = [sg, exp, genes[i], sp]
                    line = list(map(str, line))
                    print('\t'.join(line), file=fout)

    def sg_pair(self, fout=sys.stdout):
        for line in self:
            for (i1, g1), (i2, g2) in itertools.combinations(enumerate(line.genes), 2):
                if g1 is None or g2 is None:
                    continue
                id1, id2 = self.idmap[i1], self.idmap[i2]
                line = [g1, g2, '-'.join([id1, id2])]
                print('\t'.join(line), file=fout)

    def fetch(self, g1, g2):
        lines = []
        strat = False
        for line in self:
            if g1 in set(line.genes):
                strat = True
            if strat:
                lines += [line]
            if g2 in set(line.genes):
                strat = False
                return lines

    def sg_ks(self, ksfile, reflines, key='ks', max_missing=0, fout=sys.stdout, **kargs):
        d_ks = KaKsParser(ksfile).to_dict(key=key)
        for sp, lines, idx in self.parse_by_species():
            for ref, line in zip(reflines, lines):
                refg = ref.genes[0]
                if refg is None:
                    continue
                if line.missing_rate > max_missing:
                    continue
                for i, gene in enumerate(line.genes):
                    pair = (refg, gene)
                    if pair not in d_ks:
                        continue
                    ks = d_ks[pair]
                    if ks is None or ks < 0:
                        continue
                    try:
                        sg = self.idmap[idx[i]]
                    except:
                        sg = 'sg{}'.format(i)
                    line = [sg, ks, gene, sp]
                    line = list(map(str, line))
                    print('\t'.join(line), file=fout)

    def genetrees(self, pep, cds=None, outdir='genetrees', genetrees='sg.genetrees',
                  prefix=None, alignments='sg.alignments', max_missing=0.5, ncpu=50, **kargs):
        if prefix:
            genetrees = f'{prefix}.genetrees'
            alignments = f'{prefix}.alignments'
        self.trimal_opts, self.iqtree_opts = '-gappyout', '-mset GTR'
        mafft_template = '. ~/.bashrc; mafft --auto {} > {} 2> /dev/null'
        pal2nal_template = 'pal2nal.pl -output fasta {} {} > {}'
        trimal_template = 'trimal %s -in {} -out {} > /dev/null' % (
            self.trimal_opts, )
        iqtree_template = 'iqtree2 -redo -s {} %s -B 1000 -nt 1 {} > /dev/null' % (
            self.iqtree_opts, )

        mkdirs(outdir)
        self.tmpdir = outdir
        d_idmap = {}
        d_pep = seq2dict(pep)
        d_cds = seq2dict(cds)
        logger.info('Loaded {} pep, {} cds sequences'.format(
            len(d_pep), len(d_cds)))
        cdsTreefiles, cdsAlnfiles = [], []
        cmd_list = []
        for i, line in self.indexed():
            genes = [g for g in line.genes if g is not None]
            if len(genes) < 4:
                continue
            if line.missing_rate > max_missing:
                continue
            ogid = 'XG{}'.format(i)
            pepSeq = '{}/{}.pep'.format(self.tmpdir, ogid)
            cdsSeq = '{}/{}.cds'.format(self.tmpdir, ogid)
            f_pep = open(pepSeq, 'w')
            f_cds = open(cdsSeq, 'w')
            for gene, colname in zip(line.genes, self.colnames):
                if gene is None:
                    continue
                rc = d_pep[gene]
                xid = colname + '.' + gene
                rc.id = format_id_for_iqtree(xid)
                d_idmap[rc.id] = colname
                SeqIO.write(rc, f_pep, 'fasta')
                rc = d_cds[gene]
                rc.id = format_id_for_iqtree(xid)
                SeqIO.write(rc, f_cds, 'fasta')
            f_pep.close()
            f_cds.close()
            pepAln = pepSeq + '.aln'
            cdsAln = cdsSeq + '.aln'
            pepTrim = pepAln + '.trimal'
            cdsTrim = cdsAln + '.trimal'
            pepTreefile = pepTrim + '.treefile'
            cdsTreefile = cdsTrim + '.treefile'
            treefile = cdsTreefile
            cmd = '[ ! -s {} ]'.format(treefile)
            cmds = [cmd]
            cmd = mafft_template.format(pepSeq, pepAln)
            cmds += [cmd]
            iqtree_opts0 = ''  # ' -o {} '.format(root) if root else ''
            pep = True
            iqtree_opts = iqtree_opts0
            cmd = pal2nal_template.format(pepAln, cdsSeq, cdsAln)
            cmds += [cmd]
            cmd = trimal_template.format(cdsAln, cdsTrim)
            cmds += [cmd]
            cmd = iqtree_template.format(cdsTrim, iqtree_opts)
            cmds += [cmd]
            cdsTreefiles += [cdsTreefile]
            cdsAlnfiles += [cdsTrim]
            cmds = ' && '.join(cmds)
            cmd_list += [cmds]

        logger.info('Filtered {} total SGs -> {} kept (min_genes>=4, max_missing<={})'.format(
            i, len(cdsAlnfiles), max_missing))
        nbin = 10
        cmd_file = '{}/{}.cmds.list'.format(self.tmpdir, 'cds')
        run_job(cmd_file, cmd_list=cmd_list, tc_tasks=ncpu, by_bin=nbin,
                fail_exit=False)
        treefiles = cdsTreefiles
        ColinearGroups().cat_genetrees(treefiles, genetrees, idmap=d_idmap,
                                       plain=False, format_confidence='%d')
        with open(alignments, 'w') as fout:
            catAln(cdsAlnfiles, fout, idmap=d_idmap)

    def append_cols(self, seqfile, sample, prefix=None, fout=None, **kargs):
        '''Append assembled sequences as a set of new columns (one per copy).

        seqfile: HybPiper retrieve output. Sequence id = '{sample}-{colnum}_XG{row}',
        the bait id written by compmap: colnum = column number in the colnames
        file, row = alignment row number (same XG numbering as genetrees).
        New column '{sample}_{copy}' takes that sequence at row 'row' of column
        'colnum'; unmatched rows get '.'.
        prefix: output prefix -> {prefix}.alignment.csv + {prefix}.rename.list
        (default: csv to stdout, no rename.list).'''
        from Bio import SeqIO
        import re
        copy_sep = kargs.pop('copy_sep', '_')
        lines = list(self)
        pat = re.compile(r'^(\d+)_XG(\d+)$')
        d_copy = {}  # copy -> {row_index: seq id}
        nseq = 0
        for rc in SeqIO.parse(seqfile, 'fasta'):
            # id = '{sample}-{colnum}_XG{row}' (the HybPiper gene directory name)
            gname = rc.id[len(sample) + 1:] if rc.id.startswith('{}-'.format(sample)) else rc.id
            m = pat.match(gname)
            if m is None:
                logger.warning('{}: unexpected id {}, skipped'.format(seqfile, rc.id))
                continue
            col, row = int(m.group(1)), int(m.group(2))
            if not 1 <= col <= len(self.colnames) or not 1 <= row <= len(lines):
                logger.warning('{}: {} out of range, skipped'.format(seqfile, rc.id))
                continue
            if lines[row - 1].genes[col - 1] is None:
                logger.warning('{}: {} points to an empty cell, skipped'.format(seqfile, rc.id))
                continue
            copy = self.colnames[col - 1].rsplit(copy_sep, 1)[-1]
            d_copy.setdefault(copy, {})
            if row - 1 in d_copy[copy]:
                logger.warning('{}: duplicate copy {} at row {}, keep first'.format(
                    seqfile, copy, row - 1))
                continue
            d_copy[copy][row - 1] = rc.id  # full seq id, written into the new cell
            nseq += 1
        if not d_copy:
            raise ValueError('No sequences matched reference cells in {}'.format(seqfile))
        logger.info('Loaded {} seqs, copies: {}'.format(nseq,
            {c: len(g) for c, g in sorted(d_copy.items())}))
        # output: new columns {sample}_{copy} appended after the reference columns
        if prefix:
            fout = open('{}.alignment.csv'.format(prefix), 'w')
            f_colnames = open('{}.rename.list'.format(prefix), 'w')
        elif fout is None:
            fout = sys.stdout
            f_colnames = None
        else:
            f_colnames = None
        copies = sorted(d_copy)
        newcols = ['{}_{}'.format(sample, c) for c in copies]
        if f_colnames is not None:
            for i, cn in enumerate(self.colnames):
                print('{}\t{}'.format(i + 1, cn), file=f_colnames)
            for j, cn in enumerate(newcols):
                print('{}\t{}'.format(len(self.colnames) + j + 1, cn), file=f_colnames)
            f_colnames.close()
        for li, line in enumerate(lines):
            genes = list(line.genes)
            for c in copies:
                seqid = d_copy[c].get(li)
                genes.append(seqid if seqid else None)
            cells = [g if g is not None else '.' for g in genes]
            print(','.join(cells), file=fout)
        if prefix:
            fout.close()
            logger.info('Written {}.alignment.csv and {}.rename.list'.format(prefix, prefix))
    def compmap(self, cds, read1, read2=None, ref=None, prefix=None,
                ncpu=8, **kargs):
        '''Competitive mapping of reads to reference CDS baits.

        bait = all CDS of species ref (from the alignment csv), one sequence per
        gene: duplicated baits of identical sequence would compete with each
        other, which is not what the competitive mapping is meant to test.
        Reads are mapped with bwa mem and filtered to primary-only alignments
        (samtools view -G 12 -F 2304), so each read pair is assigned to at most
        one reference gene. Only the fully-retained baits are kept in the output
        bam (the all-bait bam is never written to disk) and depth is summarized
        per retained gene and copy.

        Bait ids are '{ref}-{colnum}_XG{row}', i.e. the '{taxon}-{gene}' layout
        HybPiper requires (gene id = the HybPiper gene directory, hence unique
        and '-' free), using the same XG numbering as genetrees. The bam is
        therefore ready for HybPiper reads_first.py --bwa --no-blast, which
        wants it named '{prefix}.bam' in the working directory, and the bwa
        index is dropped once the mapping is done.

        cds: all-species cds fasta (id = sp|geneid)
        read1/read2: fastq(.gz) reads
        ref: reference species name; column copies must end with _<copy>
        prefix: output prefix -> {prefix}.bait.fa {prefix}.bam {prefix}.depth.tsv
                {prefix}.target.cds'''
        from Bio import SeqIO
        import subprocess
        copy_sep = kargs.pop('copy_sep', '_')
        cov_cutoff = kargs.pop('cov_cutoff', None)
        idx_copy = []  # (col_index, copy)
        for i, cn in enumerate(self.colnames):
            if cn.startswith('{}{}'.format(ref, copy_sep)):
                idx_copy.append((i, cn[len(ref) + len(copy_sep):]))
        if not idx_copy:
            raise ValueError('No columns for ref {} in {}'.format(ref, self.colnames))
        # one cell per reference gene, preferring fully-retained rows so that a
        # gene keeps the same id as its target bait below
        d_cell = {}      # refgene -> (row, col_index), first occurrence in a retained row
        d_cell_all = {}  # refgene -> (row, col_index), first occurrence anywhere
        retained_rows = set()
        for r, line in self.indexed():
            keep = all(line.genes[i] is not None for i, _ in idx_copy)
            if keep:
                retained_rows.add(r)
            for i, _ in idx_copy:
                g = line.genes[i]
                if g is None:
                    continue
                if g not in d_cell_all:
                    d_cell_all[g] = (r, i)
                if keep and g not in d_cell:
                    d_cell[g] = (r, i)
        for g, c in d_cell_all.items():
            d_cell.setdefault(g, c)
        refs = set(d_cell_all)
        d_retain = {g: c for g, c in d_cell.items() if c[0] in retained_rows}
        logger.info('ref {} copies: {}; all ref genes: {}; fully-retained genes: {}'.format(
            ref, [c for _, c in idx_copy], len(refs), len(d_retain)))
        # ref itself may contain '-', so it is sanitized for the taxon slot; gene
        # ids never enter the id, so they may contain any character
        taxon = ref.replace('-', '_')
        def xid(cell):
            r, i = cell
            return '{}-{}_XG{}'.format(taxon, i + 1, r)
        # write bait fasta: all ref genes, one sequence per gene
        baitfile = '{}.bait.fa'.format(prefix)
        d_bait_gene = {}  # bait id (bam RNAME) -> refgene
        s_written = set()
        n = 0
        with open(baitfile, 'w') as f:
            for rc in SeqIO.parse(cds, 'fasta'):
                gene = rc.id
                if gene not in d_cell or gene in s_written:
                    continue
                bid = xid(d_cell[gene])
                rc.id = bid
                rc.description = bid
                SeqIO.write(rc, f, 'fasta')
                d_bait_gene[bid] = gene
                s_written.add(gene)
                n += 1
        logger.info('Written {} baits to {}'.format(n, baitfile))
        if n < len(refs):
            logger.warning('{} ref genes not found in {}'.format(len(refs) - n, cds))
        # bwa index + mem + primary-only + retained-only. The all-bait readings are
        # filtered on the fly, so the big all-bait bam is never written: only the
        # retained baits (the ones HybPiper assembles) reach the disk
        bamfile = '{}.bam'.format(prefix)
        idfile = '{}.kept.ids'.format(prefix)
        with open(idfile, 'w') as f:
            for g in d_retain:
                f.write('{}\n'.format(xid(d_retain[g])))
        # stale index of a previous run would not match the new bait ids
        run_cmd('rm -f {0}.amb {0}.ann {0}.bwt {0}.pac {0}.sa'.format(baitfile))
        cmds = ['bwa index {}'.format(baitfile)]
        reads = '{} {}'.format(read1, read2) if read2 else read1
        # keep bam unsorted (R1/R2 adjacent, required by HybPiper distribute); the sort
        # needed by samtools depth is piped, no sorted copy on disk
        # middle stage emits SAM text (no -b): the awk filter below reads text, a BAM
        # stream makes it emit garbage and the final samtools view -b - fails silently
        cmds.append("bwa mem -t {} {} {} | samtools view -G 12 -F 2304 -h -S - | "
                    "awk 'NR==FNR{{m[$1]=1; next}} /^@/{{print; next}} ($3 in m){{print}}' {} - | "
                    "samtools view -b - > {}".format(ncpu, baitfile, reads, idfile, bamfile))
        for c in cmds:
            logger.info(c)
            run_cmd(c, logger=logger, fail_exit=True)
        # the bwa index stores the sequence names and is of no use downstream
        run_cmd('rm -f {0}.amb {0}.ann {0}.bwt {0}.pac {0}.sa'.format(baitfile))
        os.remove(idfile)
        # depth per gene: mean depth over gene length (retained baits only, hence
        # every row of the table is retained=1)
        depthfile = '{}.depth.tsv'.format(prefix)
        run_cmd("samtools sort -@ {} {} | samtools depth -a - | awk '{{s[$1]+=$3; n[$1]++}} END {{for (g in s) print g, s[g]/n[g]}}' > {}.gene".format(
            ncpu, bamfile, depthfile))
        # summarize per copy (bam RNAME = bait id -> refgene)
        d_depth = {}
        for line in open('{}.gene'.format(depthfile)):
            g, d = line.split()
            d_depth[d_bait_gene.get(g, g)] = float(d)
        with open(depthfile, 'w') as f:
            print('gene\tcopy\tretained\tdepth', file=f)
            for g, (r, i) in sorted(d_retain.items(), key=lambda x: (x[1][1], x[0])):
                print('{}\t{}\t{}\t{}'.format(
                    g, self.colnames[i][len(ref) + len(copy_sep):], 1, d_depth.get(g, 0)), file=f)
        os.remove('{}.gene'.format(depthfile))
        # target baits for HybPiper: retained copies, ids identical to the RNAMEs
        bait_ret = '{}.target.cds'.format(prefix)
        with open(bait_ret, 'w') as f:
            for rc in SeqIO.parse(cds, 'fasta'):
                if rc.id not in d_retain:
                    continue
                bid = xid(d_retain[rc.id])
                rc.id = bid
                rc.description = bid
                SeqIO.write(rc, f, 'fasta')
        logger.info('Written {} (retained-copy baits)'.format(bait_ret))
        # summary: mean depth per copy (retained only)
        d_sum = {}
        for g, (r, i) in d_retain.items():
            d = d_depth.get(g, 0)
            d_sum.setdefault(self.colnames[i][len(ref) + len(copy_sep):], []).append(d)
        for c, ds in sorted(d_sum.items()):
            logger.info('copy {}: mean depth {:.2f} over {} genes'.format(
                c, sum(ds) / len(ds), len(ds)))
        # 每个亚基因组的深度分布：小提琴 + 中值（默认只用完全保留位点、除以总平均）
        plot_depth_violin(depthfile, outfigs='{}.depth.pdf'.format(prefix),
                          ref=ref, sample=os.path.basename(prefix.rstrip('/')))

    def win_trees(self, pep, gff, win_size=50, win_step=10, min_ratio=0.7,
                  min_seqs=4, tmpdir='tmp', exclude_ref=False, **kargs):
        mkdirs(tmpdir)
        logger.info('vars: {}'.format(vars()))
        prefix = '{}/{}'.format(tmpdir, os.path.basename(self.alignment))
        species = self.species
        logger.info('species: {}'.format(species))

        logger.info('loading {}'.format(pep))
        d_seqs = {rc.id: rc for rc in SeqIO.parse(pep, 'fasta')}
        logger.info('loading {}'.format(gff))
        d_genes = Gff(gff).get_indexed_genes()
        logger.info('loading {}'.format(self.alignment))
        lines = []
        d_alnfiles = {}
        cmd_list = []
        nolost = 0
        nsp = len(species) - 1 if exclude_ref else len(species)
        for line in self:
            line.species = species
            valid = line.nonmissing
            if len(valid) == len(species):
                nolost += 1
            if exclude_ref:
                valid = valid[1:]
            if len(valid) < min_seqs or 1.0*len(valid)/nsp < min_ratio:
                continue
            id = line.id
            outSeq = '{}.{}'.format(prefix, id)
            f = open(outSeq, 'w')
            for gene, sp in valid:
                rc = d_seqs[gene]
                rc.id = sp
                SeqIO.write(rc, f, 'fasta')
            f.close()
            alnSeq = outSeq + '.aln'
            alnTrim = alnSeq + '.trimal'
            d_alnfiles[id] = alnTrim
            cmds = []
            cmd = 'mafft --auto "{}" > "{}" 2> /dev/null'.format(
                outSeq, alnSeq)
            cmds += [cmd]
            cmd = 'trimal -automated1 -in "{}" -out "{}" &> /dev/null'.format(
                alnSeq, alnTrim)
            cmds += [cmd]
            cmds = ' && '.join(cmds)
            cmd_list += [cmds]
            lines += [line]
        logger.info('{} anchor genes'.format(len(d_alnfiles)))
        logger.info('{} non-lost genes'.format(nolost))
        if cmd_list:
            cmd_file = '{}.aln-cmds.list'.format(prefix)
            run_job(cmd_file, cmd_list=cmd_list, tc_tasks=100, by_bin=20)

        root = species[1] if exclude_ref else species[0]
        opts = '-o {}'.format(root)
        cmd_list = []
        d_treefiles = OrderedDict()
        for chrom, xlines in itertools.groupby(lines, lambda x: d_genes[x.id].chrom):
            xlines = list(xlines)
            ngenes = len(xlines)
            for i in range(0, ngenes, win_step):
                if ngenes - i < win_step:
                    continue
                id = '{}.{}'.format(chrom, i)
                lines = xlines[i:i+win_size]
                try:
                    start, end = xlines[i], xlines[i+win_step]
                except IndexError:
                    continue

                alnfiles = [d_alnfiles[line.id]
                            for line in lines if line.id in d_alnfiles]
                outALN = '{}.{}.aln'.format(prefix, id)
                with open(outALN, 'w') as f:
                    catAln(alnfiles, f)
                alnTrim = outALN + '.trimal'
                iqtreefile = alnTrim + '.treefile'
                treefile = alnTrim + '.tre'
                cmds = []
                cmd = 'trimal -automated1 -in "{}" -out "{}" &> /dev/null'.format(
                    outALN, alnTrim)
                cmds += [cmd]
                cmd = 'iqtree2 --redo -s "{}" -T 1 -B 1000 {} --mset JTT &> /dev/null'.format(
                    alnTrim, opts)
                cmds += [cmd]
                cmd = 'nw_reroot "{intre}" {root} | nw_prune - {root} '.format(
                    intre=iqtreefile, root=root,)
                cmd += ' | nw_topology -I - | nw_order - | nw_order - -c d | \
nw_order - > {}'.format(treefile)
                cmds += [cmd]
                cmds = ' && '.join(cmds)
                d_treefiles[(chrom, start, end)] = treefile
                cmd_list += [cmds]
        logger.info('{} windows'.format(len(d_treefiles)))
        cmd_file = '{}.concat-cmds.list'.format(prefix)
        run_job(cmd_file, cmd_list=cmd_list, tc_tasks=100, by_bin=1)

        outlist = '{}.wintrees'.format(self.alignment)
        outfig = outlist + '.pdf'
        from topl_distribution import Block, plot_chrom
        d_blocks = OrderedDict()
        topos = []
        f = open(outlist, 'w')
        line = ['chrom', 'istart', 'iend', 'gstart', 'gend', 'topology']
        print('\t'.join(map(str, line)), file=f)
        for (chrom, start, end), treefile in list(d_treefiles.items()):
            sid, eid = start.id, end.id
            sg, eg = d_genes[sid], d_genes[eid]
            sidx, eidx = sg.index, eg.index
            topo = ColinearGroups().get_topology(treefile)
            line = [chrom, sidx, eidx, sid, eid, topo]
            print('\t'.join(map(str, line)), file=f)
            block = Block([d_genes[sid], d_genes[eid]])
            block.topl = topo
            try:
                d_blocks[chrom] += [block]
            except KeyError:
                d_blocks[chrom] = [block]
            topos += [topo]
        f.close()

        plot_chrom(d_blocks, outfig, topos=topos)


class AlignmentLine:
    def __init__(self, line, indice=None):
        genes = []
        species = []
        for i, gene in enumerate(line.strip().split(',')):
            if i == 0:
                self.id = gene
            if indice and i not in set(indice):
                continue
            gene = None if gene == '.' or not gene else gene
            sp = None if gene is None else gene.split('|')[0]
            genes += [gene]
            species += [sp]
        self.genes = genes
        self.species = species

    def __iter__(self):
        return iter(self.genes)

    def __len__(self):
        return len(self.genes)

    def __hash__(self):
        return hash(self.id)

    @property
    def nonmissing(self):
        return [(g, sp) for g, sp in zip(self.genes, self.species) if g is not None]

    @property
    def missing_rate(self):
        return 1 - 1.0 * len(self.nonmissing) / len(self.genes)


def main():
    import sys
    subcmd = sys.argv[1]
    kargs = parse_kargs(sys.argv)
    alignment = sys.argv[2]
    if subcmd == 'wintrees':
        alignment, pep, gff = sys.argv[2:5]
        Alignment(alignment, **kargs).win_trees(pep, gff, **kargs)
    elif subcmd == 'sg_exp':
        alignment, tpmfile = sys.argv[2:4]
        Alignment(alignment, **kargs).sg_exp(tpmfile, **kargs)
    elif subcmd == 'sg_ks':
        alignment, ksfile, ref = sys.argv[2:5]
        reflines = Alignment(alignment, indice=ref)
        Alignment(alignment, **kargs).sg_ks(ksfile, reflines, **kargs)
    elif subcmd == 'sg_pair':
        alignment, idmap = sys.argv[2:4]
        Alignment(alignment, idmap=idmap, **kargs).sg_pair()
    elif subcmd == 'besthits':
        alignment, pol_indice, dip_indice, labels = sys.argv[2:6]
        blasts = sys.argv[6:]
        stats_besthits(alignment, pol_indice, dip_indice, labels,  blasts)
    elif subcmd == 'color_jcvi':
        alignment, bed,  gff = sys.argv[2:5]
        simples = sys.argv[5:]
        set_color_for_jcvi(alignment, bed,  gff, simples)
    elif subcmd == 'stat_sgs':
        alignment = sys.argv[2]
        Alignment(alignment, **kargs).stat_sgs()
    elif subcmd == 'to_bimat':  # binary matrix for iqtree
        alignment = sys.argv[2]
        Alignment(alignment, **kargs).to_bimat()
    elif subcmd == 'stat_sps':
        alignment = sys.argv[2]
        Alignment(alignment, **kargs).stat_sps()
    elif subcmd == 'high_retain':
        Alignment(alignment, **kargs).get_high_retain(**kargs)
    elif subcmd == 'high_loss':
        Alignment(alignment, **kargs).get_high_loss(**kargs)
    elif subcmd == 'genetrees':
        pep, cds = sys.argv[3:5]
        Alignment(alignment, **kargs).genetrees(pep, cds, **kargs)
    elif subcmd == 'append_cols':
        seqfile, sample = sys.argv[3:5]
        Alignment(alignment, **kargs).append_cols(seqfile, sample, **kargs)
    elif subcmd == 'compmap':
        cds, read1 = sys.argv[3:5]
        read2 = sys.argv[5] if len(sys.argv) > 5 and not sys.argv[5].count('=') else None
        Alignment(alignment, **kargs).compmap(cds, read1, read2, **kargs)
    elif subcmd == 'depthplot':
        # 只画图：soi-wgdi depthplot <深度表> [输出.pdf] [ref=物种] [sample=名字]
        plot_depth_violin(alignment, outfigs=(sys.argv[3] if len(sys.argv) > 3 else None), **kargs)
    elif subcmd == 'plot_ak':
        akfile = sys.argv[2]
        AK(akfile).plot(**kargs)
    elif subcmd == 'anc2jcvi':
        akfile = sys.argv[2]
        gff = sys.argv[3]
        AK(akfile).to_jcvi(gff=gff, **kargs)
    else:
        raise ValueError('Unknown sub command: {}'.format(subcmd))


if __name__ == '__main__':
    main()
