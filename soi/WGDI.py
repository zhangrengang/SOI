import sys
import os
import re
import itertools
import numpy as np
from collections import Counter, OrderedDict
from Bio import SeqIO
import matplotlib.pyplot as plt
from .mcscan import Gff, ColinearGroups, KaKsParser, seq2dict, Collinearity
from .OrthoFinder import catAln, format_id_for_iqtree
from lazy_property import LazyWritableProperty as lazyproperty
from .small_tools import mkdirs, flatten, test_s, test_f, parse_kargs
from .RunCmdsMP import run_cmd, run_job, logger

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
    '''Parser for WGDI ancestor.txt'''

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
        d_index = Gff(gff).get_index()
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
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = 'Arial'

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

    def plot_dotplot(self, ax=None, d_offset={}, gff=None, align='center', xy=1,  axis='x', width=0.5, label=True, gene_axis=False, fontsize=10):
        bar = plt.bar if axis == 'y' else plt.barh
        has_lab = False
        texts = []
        for sgement in self.lazy_lines(gene_axis, gff=gff):
            if sgement.chrom not in d_offset:
                continue
            offset = d_offset[sgement.chrom] + sgement.start
            lab = sgement.label if label else None
            bar(xy, len(sgement), width, offset,
                color=sgement.color, align=align, )
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
        try:
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


class Alignment:
    '''Parser for WGDI -a output (alignment file)'''

    def __init__(self, alignment, indice=None, idmap=None, colnames=None, **kargs):
        self.alignment = alignment
        self.indice = self.parse_indice(indice)  # cut columns
        if idmap is None and colnames is not None:
            idmap = colnames
        self.idmap = self.parse_idmap(idmap)  # idx - name * unique ID
        self.colnames = self.parse_colnames(colnames)  # names    * repeat allowed

    def __iter__(self):
        return self._parse()

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
                  alignments='sg.alignments'):
        self.trimal_opts, self.iqtree_opts = '-automated1', '-mset GTR'
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
        cdsTreefiles, cdsAlnfiles = [], []
        cmd_list = []
        i = 0
        for line in self:
            i += 1
            genes = [g for g in line.genes if g is not None]
            if len(genes) < 4:
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

        nbin = 10
        ncpu = 50
        cmd_file = '{}/{}.cmds.list'.format(self.tmpdir, 'cds')
        treefiles = cdsTreefiles
        ColinearGroups().cat_genetrees(treefiles, genetrees, idmap=d_idmap,
                                       plain=False, format_confidence='%d')
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
        Alignment(alignment, **kargs).genetrees(pep, cds)
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
