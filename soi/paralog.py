# coding: utf-8
"""
paralog: output paralogous gene pairs from HOGs.

For each HOG, genes from the same species form paralog pairs.
"""
import sys
from itertools import combinations
from collections import defaultdict
from .hog import HOG
from .RunCmdsMP import logger


class Paralog:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None,
				 paralog=False, nodes=None, species=None, output=None):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.paralog = paralog
		self.nodes = set(nodes) if nodes else None
		self.species = set(species) if species else None
		self.output = output

	def run(self):
		hog = HOG(ogfile=self.ogfile, orthfiles=self.orthfiles,
				  sptreefile=self.sptreefile, paralog=self.paralog)
		all_hogs = hog.pipe(write_tsv=False)
		logger.info('Loaded {} HOGs'.format(len(all_hogs)))

		fout = open(self.output, 'w') if self.output else sys.stdout
		count = 0
		for hog_id, hog_rec in all_hogs.items():
			node_id = hog_rec['node_id']
			if self.nodes and node_id not in self.nodes:
				continue
			genes = hog_rec['genes']
			sp_genes = defaultdict(list)
			for g in genes:
				sp = g.split('|')[0] if '|' in g else g.rsplit('_', 1)[0]
				sp_genes[sp].append(g)
			for sp, gs in sp_genes.items():
				if self.species and sp not in self.species:
					continue
				if len(gs) < 2:
					continue
				for g1, g2 in combinations(sorted(gs), 2):
					fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
						g1, g2, node_id, sp, hog_id))
					count += 1

		if self.output:
			fout.close()
		logger.info('Output {} paralog pairs'.format(count))


def xmain(**kargs):
	Paralog(**kargs).run()
