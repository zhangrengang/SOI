#coding: utf-8
import sys, os
import argparse
from math import sqrt, ceil
import networkx as nx
import numpy as np
import itertools
import matplotlib.pyplot as plt
from mcscan import Collinearity, Gff

__version__ = '0.1'
__LastModified__ = '20200616'
__Example__ = None
def add_ploidy_opts(parser):
	parser.add_argument('--window_size', type=int, default=50, help="window_size. default=%(default)s")
	parser.add_argument('--window_step', type=int, default=10, help="window_step. default=%(default)s")
	parser.add_argument('--min_block', type=int, default=10, help="min genes for a block. default=%(default)s")
	parser.add_argument('--max_distance', type=int, default=20, help="max_distance. default=%(default)s")
	parser.add_argument('--max_ploidy',  type=int, default=10, help="x limit. default=%(default)s")
	parser.add_argument('--min_overlap', type=float, default=20, help="min_overlap. default=%(default)s")
	parser.add_argument('--color',  type=str, default=None, help="bar fill color. default=%(default)s")
	parser.add_argument('--edgecolor', type=str, default=None, help="bar edge color. default=%(default)s")
	
def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter,\
		epilog="Version: {}\nLast Modification Date: {}".format(__version__,__LastModified__),\
		version="Version: {}".format(__version__),\
		description="Example: {}".format(__Example__))
	parser.add_argument('-s', '--collinearity', metavar='INPUT_BLOCK_FILE', type=str, required=True, help="the blocks (*.collinearity, output of MCSCANX)")
	parser.add_argument('-g', '--gff', metavar='INPUT_GENE_GFF_FILE', type=str, required=True, help="the annotation gff file (one of MCSCANX input)")
	parser.add_argument('-r', '--ref', metavar='reference', type=str, required=True, help="reference species")
	parser.add_argument('-q', '--qry', metavar='queries', nargs='+', type=str, required=True, help="query species")
	parser.add_argument('-o', '--output', metavar='OUTPUT_FILE_PREFIX', type=str, default=None, help="the output file prefix.")
	parser.add_argument('--format', metavar='figure file out format', action='append', default=['pdf', 'png'], help="default=%(default)s")
	parser.add_argument('--nrow', metavar='nrow', type=int, default=None, help="number of rows. default=%(default)s")
	parser.add_argument('--min_same_block', type=int, default=25, help="min gene number in a block on the same chromosome. default=%(default)s")
	add_ploidy_opts(parser)
	args = parser.parse_args()
	sps = [args.ref] + args.qry
	if args.output is None:
		sps = map(lambda x:x[:2], sps)
		args.output = '-'.join(sps) + '_' + str(args.window_size)
	if args.nrow is None:
		args.nrow = int(ceil(sqrt(len(args.qry))))
	if args.window_step is None:
		args.window_step = args.window_size / 5
	if args.min_overlap is None:
		args.min_overlap = args.window_size / 2.5
	else:
		args.min_overlap = args.min_overlap*args.window_size
	args.ncol = int(ceil(1e0*len(args.qry) / args.nrow))
	args.outfigs = [args.output+'.'+fmt for fmt in args.format]
	suptitle = 'Reference: ' + args.ref
	xlabel = 'Relative Ploidy by {}-gene windows'.format(args.window_size)
	args.suptitle = '{} ({})'.format(xlabel, suptitle)
	args.titles = args.qry
	print >>sys.stderr, '{} x {} figure'.format(args.nrow, args.ncol)
	return args
def main():
	args = makeArgparse()
	plot_fold(**args.__dict__)
def plot_fold(collinearity, gff, ref, qry, **kargs):
	d_ortholog_graph, d_paralog_graph = parse_collinearity(collinearity, ref, qry, **kargs)
	d_coord_path, d_coord_graph = parse_gff(gff, [ref]+qry)
	data = []
	for sp in qry:
		d_fold = get_ploidy(d_coord_path[ref], d_coord_graph[ref],
					d_coord_graph[sp], d_ortholog_graph[sp], 
					d_paralog_graph[sp], **kargs)
		data += [np.array(sorted(d_fold.items()))]
	return plot_bars(data, **kargs)
def plot_bars(data, titles, ax=None, outfigs=None, nrow=1, ncol=1, suptitle=None, max_ploidy=10, 
			color='white', edgecolor='black', ylabel='Number of windows', 
			xlabel='Relative ploidy per reference block', **kargs):
	if ax is None:
		if nrow*ncol == 1:
			ax = plt.subplot(111)
			ax = [ax]
		else:
			fig, ax = plt.subplots(nrow, ncol, sharex=True, figsize=(10*ncol/2, 8*ncol/2))
			cells = list(itertools.product(range(nrow), range(ncol)))
			ax = [ax[cell] for cell in cells]
	else:
		ax = [ax]
	tick_label = range(0, max_ploidy+1)
	for i, (dat, title, sax) in enumerate(zip(data, titles, ax)):
		try:
			x = dat[:, 0]
			y = dat[:, 1]
			sax.bar(x, y, align='center', color=color, edgecolor=edgecolor)
		except IndexError:
			pass
		if title is not None:
			sax.set_title(title)
		sax.set_xlim(0, max_ploidy)
		if xlabel is not None and i >= (nrow-1)*ncol:
			sax.set_xlabel(xlabel)
		if ylabel is not None and i%ncol == 0:
			sax.set_ylabel(ylabel)
	plt.xticks(tick_label)
	if suptitle is not None:
		plt.suptitle(suptitle)
	#plt.subplots_adjust()
	#plt.xlabel(xlabel)
	if outfigs is not None:
		for outfig in outfigs:
			plt.savefig(outfig)
	else:
		return ax
def parse_collinearity(collinearity, ref, qry, min_block=10, min_same_block=25, **kargs):
	d_ortholog_graph = {}
	d_paralog_graph = {}
	#d_ref_ortholog = {}
	for sp in qry:
		d_ortholog_graph[sp] = nx.Graph()
		d_paralog_graph[sp] = nx.Graph()
	#	d_ref_ortholog[sp] = {}
	
	qry = set(qry)
	for rc in Collinearity(collinearity):
		if rc.chr1 == rc.chr2 and rc.N < min_same_block:
			continue
		if rc.N < min_block:
			continue
		sp1, sp2 = rc.species
		if sp1 == sp2 and sp1 in qry:
			d_paralog_graph[sp1].add_edges_from(rc.pairs)
			continue
		elif (sp1 == ref and sp2 in qry):
	#		r,q = sp1, sp2
			d_ortholog_graph[sp2].add_edges_from(rc.pairs)
		elif (sp2 == ref and sp1 in qry):
	#		r, q = sp2, sp1
			d_ortholog_graph[sp1].add_edges_from(rc.pairs)
	#	else:
	#		continue
	#	try: d_ref_ortholog[sp][r] += [q]
	#	except KeyError: d_ref_ortholog[sp][r] = [q]
	# ortholog -> paralog
	#for sp in qry:
	#	for r, qs in d_ref_ortholog[sp].items():
	#		for q1, q2 in itertools.combinations(qs, 2):
	#			d_paralog_graph[sp].add_edge(q1, q2)
	return d_ortholog_graph, d_paralog_graph
	
def parse_gff(gff, sps):
	d_coord_graph = {}
	for sp in sps:
		d_coord_graph[sp] = nx.Graph()
	sps = set(sps)
	d_gff = {}
	for line in Gff(gff):
		if not line.species in sps:
			continue
		key = (line.species, line.chrom)
		try: d_gff[key] += [line]
		except KeyError: d_gff[key] = [line]

	d_coord_path = {}
	for (sp, chrom), lines in d_gff.items():
		lines = sorted(lines, key=lambda x:(x.start, -x.end))
		genes = [line.gene for line in lines]
		try: d_coord_path[sp] += [genes]
		except KeyError: d_coord_path[sp] = [genes]
		for i, line in enumerate(lines):
			d_coord_graph[sp].add_node(line.gene, chrom=chrom, index=i)
		for i, line in enumerate(lines[1:]):
			edge = (lines[i].gene, line.gene)
			d_coord_graph[sp].add_edge(*edge)
	return d_coord_path, d_coord_graph
def get_ploidy(ref_coord_paths, ref_coord_graph, qry_coord_graph, rq_ortholog_graph, 
			qry_paralog_graph=None, # not used
			window_size=20, window_step=10, **kargs):
	'''for each query. 参考分段，每段对应query有多少块'''
	
	d_fold = {}
	for path in ref_coord_paths:
		for i in range(0, len(path), window_step):
			start, end = i, i+window_size
			if end > len(path):
				end = len(path)
			if end - start < window_size/2:
				continue
			bin = path[start:end]
			orthologs = []
			for gene in bin:
				if rq_ortholog_graph.has_node(gene):
					orthologs += rq_ortholog_graph.neighbors(gene)
			if len(orthologs) <2:	# ploidy=0
				continue
			qry_clusters = cluster_genes(orthologs, qry_coord_graph)
			qry_blocks = list(nx.connected_components(qry_clusters))
			if not qry_blocks:	# not into blocks
				continue
			ref_blocks = map_graph(bin, rq_ortholog_graph, qry_blocks)
			ref_clusters = overlap_blocks(ref_blocks, ref_coord_graph, **kargs)
			ncmpt1 = len(qry_blocks)
		#	print orthologs, qry_blocks, ref_blocks
		#	if len(ref_blocks) == 1:
		#		ncmpt2 = 1
		#	else:
			ncmpt2 = max(map(lambda x:len(x), nx.connected_components(ref_clusters)))
		#	insect_nodes = set(clusters.nodes()) & set(qry_paralog_graph.nodes())
		#	for n1, n2 in itertools.combinations(insect_nodes, 2):
		#		if qry_paralog_graph.has_edge(n1, n2):
		#			clusters.add_edge(n1, n2)
		#	ncmpt2 = len(list(nx.connected_components(clusters)))
		#	if ncmpt2 == 0:
		#		continue
		#	ploidy = int(round(1e0*ncmpt1/ncmpt2, 0)) if qry_paralog_graph else ncmpt1
		#	#ploidy = ncmpt1
			ploidy = ncmpt2
#			print ncmpt1, ncmpt2, ploidy
			try: d_fold[ploidy] += 1
			except KeyError: d_fold[ploidy] = 1
	return d_fold
def map_graph(bin, rq_ortholog_graph, qry_blocks):
	'''将qry的block映射到ref的block'''
	ref_blocks = []
	for block in qry_blocks:
		ref_block = []
		for gene in block:
			ref_block += list(set(rq_ortholog_graph[gene]) & set(bin))
		ref_blocks += [ref_block]
	return ref_blocks
def overlap_blocks(blocks, coord_graph, min_overlap=3, **kargs):
	'''将有overlap的block连起来'''
	blocks = map(tuple, blocks)
	G = nx.Graph()
	for b in blocks:
		G.add_node(b)
	for b1, b2 in itertools.combinations(blocks, 2):
		i1 = map(lambda x: coord_graph.node[x]['index'], b1)
		i2 = map(lambda x: coord_graph.node[x]['index'], b2)
		min_i1, max_i1 = min(i1), max(i1)
		min_i2, max_i2 = min(i2), max(i2)
		if min(max_i1, max_i2) - max(min_i1, min_i2) + 1 >=min_overlap: # overlap
			G.add_edge(b1, b2)
	return G
def cluster_genes(genes, coord_graph, max_distance=25, **kargs):
	'''根据坐标，将gene聚类为block'''
	d_bin = {}
	for gene in genes:
		try: chrom = coord_graph.node[gene]['chrom']
		except KeyError: continue
		try: d_bin[chrom] += [gene]
		except KeyError: d_bin[chrom] = [gene]
	G = nx.Graph()
	for chrom, genes in d_bin.items():
		genes = sorted(genes, key=lambda x: coord_graph.node[x]['index'])
		for i, gene in enumerate(genes[1:]):
			n1, n2 = genes[i], gene
			i1, i2 = coord_graph.node[n1]['index'], coord_graph.node[n2]['index']
			if i2 - i1 < max_distance:
				G.add_edge(n1, n2)
	return G

if __name__ == '__main__':
	main()
