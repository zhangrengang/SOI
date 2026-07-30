# coding: utf-8
"""
paralog: output paralogous gene pairs produced at each branch from HOGs.

At an internal node, a parent HOG splits into child HOGs.  Genes of the same
species in different child HOGs are paralogs produced on that branch.
"""
import sys
from .hog import HOG, iter_branch_paralogs
from .RunCmdsMP import logger


class Paralog:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None,
				 paralog=False, nodes=None, species=None, output=None):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.paralog = paralog
		self.nodes = nodes
		self.species = species
		self.output = output

	def run(self):
		hog = HOG(ogfile=self.ogfile, orthfiles=self.orthfiles,
				  sptreefile=self.sptreefile, paralog=self.paralog)
		all_hogs = hog.pipe(write_tsv=False)
		logger.info('Loaded {} HOGs'.format(len(all_hogs)))

		fout = open(self.output, 'w') if self.output else sys.stdout
		fout.write('#gene1\tgene2\tnode\tspecies\tHOG_id\n')
		count = 0
		for g1, g2, node_id, sp, hog_id in iter_branch_paralogs(
				all_hogs, hog.tree, self.nodes, self.species):
			fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
				g1, g2, node_id, sp, hog_id))
			count += 1

		if self.output:
			fout.close()
		logger.info('Output {} paralog pairs'.format(count))


def xmain(**kargs):
	Paralog(**kargs).run()
