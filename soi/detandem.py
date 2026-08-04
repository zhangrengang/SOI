# coding: utf-8
"""
Remove tandem duplicate genes from orthogroups.

A tandem duplicate is defined as two or more genes in the same OG that:
  1. belong to the same species and chromosome
  2. have gene index distance < tandem_dist (default 200)

When ortholog/collinearity files are provided, the gene with the highest
degree in the ortholog graph is retained. Otherwise, one gene is kept at
random from each tandem cluster.
"""

import sys
import random
from collections import defaultdict

from .OrthoFinder import OrthoMCLGroup
from .mcscan import XGff, ColinearGroups
from .RunCmdsMP import logger


class Detandem:
	"""Detect and remove tandem duplicate genes from orthogroups."""

	def __init__(self, ogfile=None, gfffile=None, orthfiles=None,
				 tandem_dist=200, out_tandem=None, **kargs):
		self.ogfile = ogfile
		self.gfffile = gfffile
		self.orthfiles = orthfiles
		self.tandem_dist = tandem_dist
		self.out_tandem = out_tandem
		self._tandem_clusters = []  # [(retained_gene, [all_genes_in_cluster]), ...]

	def run(self):
		# --- 1. parse GFF, build gene -> (chrom, index) map ---
		logger.info(f'Parsing GFF from {self.gfffile}')
		gff_parser = XGff(self.gfffile)
		self.d_genes = gff_parser.get_indexed_genes()
		logger.info(f'  {len(self.d_genes)} genes indexed')

		# --- 2. parse ortholog/collinearity files, build graph ---
		self.graph = None
		if self.orthfiles:
			logger.info(f'Loading ortholog/collinearity pairs from {self.orthfiles}')
			self.graph = ColinearGroups(self.orthfiles,
										noparalog=False).graph
			logger.info(f'  {self.graph.number_of_nodes()} genes, '
						f'{self.graph.number_of_edges()} edges')
		else:
			logger.info('No ortholog/collinearity files provided; '
						'will keep random gene per tandem cluster')

		# --- 3. iterate OGs, detect & remove tandem duplicates ---
		logger.info(f'Processing OGs from {self.ogfile}')
		total_removed = 0
		total_ogs = 0
		for og in OrthoMCLGroup(self.ogfile):
			total_ogs += 1
			filtered_genes = self._filter_tandem(og)
			removed = len(og.genes) - len(filtered_genes)
			total_removed += removed
			if filtered_genes:
				sys.stdout.write(f'{og.ogid}: {" ".join(sorted(filtered_genes))}\n')

		logger.info(f'Done. {total_removed} tandem duplicates removed from {total_ogs} OGs.')

		# --- 4. output tandem clusters if requested ---
		if self.out_tandem:
			self._write_tandem_clusters()

	def _filter_tandem(self, og):
		"""Remove tandem duplicates from one OG. Return list of retained genes."""
		spchr_genes = defaultdict(list)
		missing = []
		for gene_id in og.genes:
			gff_line = self.d_genes.get(gene_id)
			if gff_line is None:
				missing.append(gene_id)
				continue
			key = (gff_line.species, gff_line.chrom)
			spchr_genes[key].append(gff_line)

		retained = list(missing)

		for (sp, chrom), genes in spchr_genes.items():
			genes.sort(key=lambda g: g.index)
			clusters = []
			current_cluster = [genes[0]]
			for i in range(1, len(genes)):
				if genes[i].index - genes[i-1].index < self.tandem_dist:
					current_cluster.append(genes[i])
				else:
					clusters.append(current_cluster)
					current_cluster = [genes[i]]
			clusters.append(current_cluster)

			for cluster in clusters:
				if len(cluster) == 1:
					retained.append(cluster[0].gene)
				else:
					best = self._pick_best(cluster)
					retained.append(best)
					if self.out_tandem is not None:
						self._tandem_clusters.append(
							(best, [g.gene for g in cluster]))

		return retained

	def _pick_best(self, cluster):
		"""Pick the gene with highest degree in ortholog graph, or random."""
		if self.graph is None:
			return random.choice(cluster).gene
		best = max(cluster, key=lambda g: self.graph.degree(g.gene))
		return best.gene

	def _write_tandem_clusters(self):
		"""Write tandem clusters to file, retained gene marked with *."""
		with open(self.out_tandem, 'w') as fout:
			for retained, genes in self._tandem_clusters:
				line = []
				for g in genes:
					if g == retained:
						line.append(f'*{g}')
					else:
						line.append(g)
				fout.write('\t'.join(line) + '\n')
		logger.info(f'Tandem clusters written to {self.out_tandem}')
