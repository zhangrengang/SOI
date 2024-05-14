import sys
import itertools
from collections import Counter
import networkx as nx
from .mcscan import ColinearGroups, Gff, GffGraph, SyntenyGraph
from .RunCmdsMP import logger

def akr(collinearity, gff, anc, spsd, rounds=3):
	synG = ColinearGroups(collinearity, spsd=spsd, nosamechr=True, noparalog=False).to_graph()
	gffG = Gff(gff).to_graph()
	# unify node object
	synG = convert_synG(synG, gffG)
	# clean non-syntenic nodes
	nonSynNodes = set(gffG)-set(synG)
	gffG.remove_internals(nonSynNodes)
	logger.info('{} non-syntenic nodes removed'.format(len(nonSynNodes)))
	ancG = gffG.subgraph([node for node in gffG.nodes if node.species == anc])
	logger.info('{} anc nodes'.format(len(ancG)))

	TANDEM = identify_tandem(synG)
	gffG.remove_internals(TANDEM)
	logger.info('{} tandem nodes removed'.format(len(TANDEM)))

	gffG.index()

#	with open('gff.{}.gfa'.format(2), 'w') as fout:
#		gffG.to_gfa(fout)
	ancG = gffG.subgraph([node for node in gffG.nodes if node.species == anc])
	ancG = GffGraph(ancG)
	logger.info('{} anc nodes'.format(len(ancG)))

	for i in range(rounds):
		logger.info('round {}'.format(i))
#		logger.info(len(TANDEM))
		d_ancs = process_synG(synG, ancG)
#		logger.info(len(TANDEM))
		logger.info(len(ancG))
		# remove tandem
		logger.info('{} tandems'.format(len(set(ancG.nodes) - set(d_ancs))))
		print(set(ancG.nodes) - set(d_ancs))
		ancG.remove_internals(set(ancG.nodes) - set(d_ancs))
		logger.info((len(ancG), len(list(ancG.starts))))
#		with open('sny.{}.gfa'.format('x'), 'w') as fout:
#			ancG.to_gfa(fout)
		np, nn = insert_syn(ancG, gffG, d_ancs, synG, TANDEM=TANDEM)
		logger.info('inserted {} paths, {} nodes'.format(np, nn))
		logger.info((len(ancG), len(list(ancG.starts))))
#		logger.info((list(synG.nodes())[:10], list(gffG.nodes())[:10]))

		with open('sny.{}.gfa'.format(i), 'w') as fout:
			ancG.to_gfa(fout)
#	d_ancs = process_synG(synG, ancG, TANDEM=TANDEM)
#	ancG.remove_internals(set(ancG.nodes) - set(d_ancs))
	logger.info(len(ancG))
	ancG.to_wgdi(anc+'.akr')
	ancG.to_idmap()
def insert_syn(ancG, gffG, d_ancs, synG, max_dist=5, TANDEM=set([])):
	starts = list(ancG.starts)
	i,j = 0,0
	insert_paths, insert_nodes = [], []
	for start_node in starts:
		chrom = list(ancG.iter_chrom(start_node))
		for i in range(len(chrom)-1):
			start, end = chrom[i:i+2]
			sgs, egs = d_ancs[start], d_ancs[end]
			d_sgs = {sp: list(g) for sp, g in itertools.groupby(sgs, key=lambda x:x.species)}
			d_egs = {sp: list(g) for sp, g in itertools.groupby(egs, key=lambda x:x.species)}
			shared_species = set(d_sgs) & set(d_egs)
			paths = []
			for sp in shared_species:
				sg, eg = d_sgs[sp], d_egs[sp]
				for g1, g2 in itertools.product(sg, eg):
					if g1 in TANDEM or g2 in TANDEM or ({g1, g2} & {start, end}):
						continue
					if is_adj(g1, g2, max_dist, min_dist=2):
						#print(g1, g2)
						path = gffG.lazy_fetch_chrom(g1, g2)
						path.species = sp
						path.score = synG.score_path(path)
						print(i, path.score, len(path), sp, g1, g2, start, end, path)
						paths += [path]
			if not paths:
				continue
			path = max(paths, key=lambda x:x.score)
			path = path[1:-1]
			insert_paths += [(start, end, path)]
			insert_nodes += path
	d_count = Counter(insert_nodes)
	for start, end, path in insert_paths:
		path = [n for n in path if d_count[n] == 1 and n not in ancG]
		if len(path) == 0:
			continue
		ancG.insert_path(start, end, path)
		i += 1
		j += len(path)
	return i,j
def convert_synG(synG, gffG):
	newG = SyntenyGraph()
	d_gffG = {n.id:n for n in gffG.nodes}
	i,j = 0,0
	for n1, n2 in synG.edges:
		i += 1
		if not n1 in d_gffG or not n2 in d_gffG:
			continue
		j += 1
		attr = synG[n1][n2]
		n1, n2 = d_gffG[n1], d_gffG[n2]
		newG.add_edge(n1, n2, **attr)
	assert i ==j #print (i, j)
	return newG

def is_adj(g1, g2, max_dist=3, min_dist=0):
	dist = abs(g1.index - g2.index)
	if g1.species == g2.species and g1.chrom == g2.chrom and min_dist <= dist <= max_dist:
		return True
	return False
def identify_tandem(synG, max_dist=1):
	tandems = set([])
	for cmpt in nx.connected_components(synG):
		tanG = nx.Graph()
		for g1, g2 in itertools.combinations(cmpt, 2):
			if is_adj(g1, g2, max_dist):    # tandem
				tanG.add_edge(g1, g2)
		grps = list([list(c) for c in nx.connected_components(tanG)])
		for grp in grps:
			d_weight = {node: synG.score_node(node) for node in grp }
			core = max(grp, key=lambda x: d_weight[x])
			tandem = set(grp) - set([core])
			tandems = tandems | tandem
	return tandems

def process_synG(synG, ancG, max_dist=1, min_sps=2):
	d_ancs = {}
	i,j,k,m,n = 0,0,0,0,0
	for cmpt in nx.connected_components(synG):
		n +=1
		cmpt = sorted(cmpt)
		ancNodes = [nd for nd in cmpt if nd in ancG]
		sps = {nd.species for nd in cmpt}
		if len(ancNodes) ==0 or len(sps) < min_sps:
			i+=1
			continue
		elif len(ancNodes) == 1:	# only one core
			j +=1
			d_ancs[ancNodes[0]] = cmpt
			continue
		sg = synG.subgraph(cmpt)
		# multi core, to split
		m +=1
		#bin
		cores = set(ancNodes)
		d_groups = {}
		for node in cmpt:
			if node in cores:
				core = node
			else:
				dists = {core: nx.shortest_path_length(sg, node, core, weight='weight') for core in cores}
				#print(dists)
				core = min(cores, key=lambda x: dists[x])
			try: d_groups[core] += [node]
			except KeyError: d_groups[core] = [node]
		d_ancs.update(d_groups)	# >=1 groups
	return d_ancs

def test_single(synG, anc):
	for cmpt in nx.connected_components(synG):
		ancNodes = [nd for nd in cmpt if nd.species == anc]
		if len(ancNodes) == 0:
			print(cmpt)
def print_edge(synG, node):
	for xn in synG.nodes:
		if xn.id == node:
			break
	print(node, synG[xn])
def main():
	collinearity, gff, anc, spsd = sys.argv[1:5]
	akr(collinearity, gff, anc, spsd)
