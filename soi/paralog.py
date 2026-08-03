# coding: utf-8
"""
paralog: output paralogous gene pairs produced at each branch from HOGs,
         and index collinearity blocks by branch-paralog content.

Internal nodes: sibling child HOGs (same parent) -> genes of the same species
  across different child HOGs are paralogs on that branch.
Leaf species: genes of the same species within the same HOG -> paralogs.
"""
from collections import defaultdict, Counter
from .hog import HOG
from .mcscan import XCollinearity, Gff
from .RunCmdsMP import logger


# ---------------------------------------------------------------------------
#  paralog output (used directly or by indexer)
# ---------------------------------------------------------------------------

class Paralog:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None,
				 nodes=None, species=None, prefix='paralog', **kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.kargs = kargs

	def run(self):
		"""Load HOGs, write paralog TSV, return {branch: frozenset(pairs)}, root_name."""
		hog = HOG(ogfile=self.ogfile, orthfiles=self.orthfiles,
				  sptreefile=self.sptreefile, **self.kargs)
		hog.pipe(write_tsv=False)
		logger.info('Loaded {} HOGs'.format(len(hog.all_hogs)))

		fpath = self.prefix + '.paralog.tsv'
		branch_pairs, count = self._write_and_group(hog, fpath)
		logger.info('Output {} paralog pairs to {}'.format(count, fpath))
		return {b: frozenset(ps) for b, ps in branch_pairs.items()}, hog.tree.name

	def _write_and_group(self, hog, fpath):
		"""Write paralog TSV, return {branch: set of canonical (g1,g2)}, count."""
		branch_pairs = defaultdict(set)
		with open(fpath, 'w') as fout:
			fout.write('#gene1\tgene2\tnode\tspecies\tHOG_id\n')
			count = 0
			for g1, g2, node_id, sp, hog_id in hog.iter_branch_paralogs(
					self.nodes, self.species):
				fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
					g1, g2, node_id, sp, hog_id))
				pair = (g1, g2) if g1 < g2 else (g2, g1)
				branch_pairs[node_id].add(pair)
				count += 1
		return branch_pairs, count


# ---------------------------------------------------------------------------
#  branch-paralog indexer (default mode)
# ---------------------------------------------------------------------------

class ParalogIndexer:
	"""Assign collinearity blocks to branches based on paralog content."""

	def __init__(self, ogfile, orthfiles, sptreefile,
				 self_synteny=None, min_n=0, gff=None, min_dist=None,
				 pi_cutoff=0.05, nodes=None, species=None,
				 prefix='paralog_index', **hog_kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.self_synteny = self_synteny or orthfiles
		self.min_n = min_n
		self.gff = gff
		self.min_dist = min_dist
		self.pi_cutoff = pi_cutoff
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.hog_kargs = hog_kargs

		# lazy
		self._branch_pairs = None
		self._root_branch = None

	def _load_branch_pairs(self):
		"""Load HOGs and paralog pairs via Paralog."""
		if self._branch_pairs is not None:
			return
		paralog = Paralog(ogfile=self.ogfile, orthfiles=self.orthfiles,
						  sptreefile=self.sptreefile, prefix=self.prefix,
						  nodes=self.nodes, species=self.species,
						  **self.hog_kargs)
		self._branch_pairs, self._root_branch = paralog.run()

	def _compute_pi(self, block_pairs, branch_pairs):
		"""Compute Paralogue Index = |intersection| / |block_pairs|."""
		if not block_pairs:
			return 0.0
		common = sum(1 for p in block_pairs if p in branch_pairs)
		return common / len(block_pairs)

	def _canonical_pair(self, g1, g2):
		return (g1, g2) if g1 < g2 else (g2, g1)

	def _filter_tandem(self, blocks):
		"""Remove tandem blocks if -d is set.  Lazy: load GFF on first call."""
		if not self.min_dist or not self.gff:
			yield from blocks
			return
		d_gene = Gff(self.gff).get_indexed_genes()
		for rc in blocks:
			if _is_tandem(rc, d_gene, self.min_dist):
				continue
			yield rc

	def assign(self):
		"""Assign blocks to branches, return {branch: [(block, sp1, sp2, N, pi), ...]}."""
		self._load_branch_pairs()
		assigned = defaultdict(list)
		threshold = self.pi_cutoff
		root = self._root_branch

		for rc in self._filter_tandem(XCollinearity(self.self_synteny)):
			if rc.N < self.min_n:
				continue
			if rc.species1 != rc.species2:
				continue

			block_pairs = [self._canonical_pair(g1, g2) for g1, g2 in rc.pairs]

			best_branch = root
			best_pi = 0.0
			for branch, bp_set in self._branch_pairs.items():
				pi = self._compute_pi(block_pairs, bp_set)
				if pi > best_pi:
					best_pi = pi
					best_branch = branch

			if best_pi < threshold:
				best_branch = root

			assigned[best_branch].append(
				(rc.block, rc.species1, rc.species2, rc.N, best_pi))

		logger.info('Assigned blocks to {} branches'.format(len(assigned)))
		return assigned

	def write_stats(self, assigned):
		"""Write per-branch per-species statistics TSV."""
		fpath = self.prefix + '.index.stats.tsv'
		stats = defaultdict(lambda: [0, 0, 0, 0.0])
		for branch, items in assigned.items():
			for block, sp1, sp2, N, pi in items:
				for sp in (sp1, sp2):
					if sp is None:
						continue
					s = stats[(branch, sp)]
					s[0] += 1
					s[1] += N
					s[2] += int(pi * N)
					s[3] += pi

		with open(fpath, 'w') as fout:
			fout.write('#branch\tspecies\tblocks\tgene_pairs\tparalog_pairs\tmean_PI\n')
			for (branch, sp), (blocks, gp, pp, sum_pi) in sorted(stats.items()):
				mean_pi = sum_pi / blocks if blocks else 0.0
				fout.write('{}\t{}\t{}\t{}\t{}\t{:.4f}\n'.format(
					branch, sp, blocks, gp, pp, mean_pi))
		logger.info('Stats written to {}'.format(fpath))

	def write_blocks(self, assigned):
		"""Write assigned blocks per branch."""
		for branch, items in assigned.items():
			fpath = '{}.index.{}.blocks'.format(self.prefix, branch)
			with open(fpath, 'w') as fout:
				for block, sp1, sp2, N, pi in items:
					fout.write(block)
		logger.info('Block files written for {} branches'.format(len(assigned)))


def _is_tandem(rc, d_gene, min_dist):
	"""Check if block is tandem on the same chromosome within min_dist."""
	if rc.chr1 != rc.chr2:
		return False
	positions = []
	for g1, g2 in rc.pairs:
		p1 = d_gene.get(g1, {}).get('pos')
		p2 = d_gene.get(g2, {}).get('pos')
		if p1 is not None and p2 is not None:
			positions.append(abs(p1 - p2))
	if positions and max(positions) < min_dist:
		return True
	return False


# ---------------------------------------------------------------------------
#  CLI entry
# ---------------------------------------------------------------------------

def xmain(**kargs):
	no_index = kargs.pop('no_index', False)
	if no_index:
		kargs.pop('output', None)
		Paralog(**kargs).run()
	else:
		indexer = ParalogIndexer(**kargs)
		assigned = indexer.assign()
		indexer.write_stats(assigned)
		indexer.write_blocks(assigned)
