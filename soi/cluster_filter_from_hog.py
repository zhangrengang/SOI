#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从HOG信息创建单拷贝OG文件
根据HOG（Hierarchical Orthologous Groups）信息，对OG文件进行处理，
确保每个物种在每个OG中只保留一个基因。

优化版：针对大OG文件和大物种树进行了性能优化。
"""

import time
import argparse
from collections import defaultdict
from contextlib import nullcontext

# 外部依赖
import networkx as nx
from ete3 import Tree

# 假设这些是项目内部模块，请根据实际情况调整导入路径
from .hog import HOG
from .OrthoFinder import OrthoMCLGroup, gene_format_o
from .RunCmdsMP import logger
from .mcscan import ColinearGroups


def load_hog_info_from_class(ogfile, orthfiles, sptreefile, paralog=False):
	"""
	从HOG类加载HOG信息，并预计算每个HOG的物种集合与物种数（使用独立字典）。

	:param ogfile: OG文件路径
	:param orthfiles: 正交文件路径列表
	:param sptreefile: 物种树文件路径
	:param paralog: 是否包含旁系同源，默认False
	:return: (all_hogs, tree, hog_species_count)
	"""
	hog_instance = HOG(ogfile=ogfile, orthfiles=orthfiles,
					sptreefile=sptreefile, paralog=paralog)
	all_hogs = hog_instance.pipe(write_tsv=False)

	# 预计算每个HOG的物种数，存储为独立字典
	hog_species_count = {}
	for hog_id, rec in all_hogs.items():
		species_set = {gene_format_o(g)[0] for g in rec.genes}
		hog_species_count[hog_id] = len(species_set)

	return all_hogs, hog_instance.tree, hog_species_count

def build_deletion_map(all_hogs, tree, hog_species_count, gene_degree=None):
	"""
	遍历树一次，构建要删除的基因集合（优化版：标记删除 + 预排序 + 直接使用物种数）

	:param all_hogs: HOG字典
	:param tree: 物种树对象
	:param hog_species_count: 字典 {hog_id: species_count}
	:param gene_degree: 基因度数（可选，用于叶子节点基因保留策略）
	:return: 待删除的基因集合
	"""
	# 预构建节点到HOG的映射（按hog_id排序）
	node_to_hogs = defaultdict(list)
	for hog_id, hog_record in all_hogs.items():
		node_to_hogs[hog_record.node_id].append((hog_id, hog_record))
	for lst in node_to_hogs.values():
		lst.sort(key=lambda x: x[0])   # 按hog_id排序

	deleted_hog_ids = set()   # 标记删除的HOG
	deleted_genes = set()	 # 待删除的基因

	# 递归删除当前HOG及其所有后代（用于子节点比较）
	def mark_hog_and_descendants(hid):
		if hid in deleted_hog_ids:
			return
		stack = [hid]
		while stack:
			cur_id = stack.pop()
			if cur_id in deleted_hog_ids:
				continue
			cur_rec = all_hogs.get(cur_id)
			if not cur_rec:
				continue
			deleted_hog_ids.add(cur_id)
			deleted_genes.update(cur_rec.genes)
			for child_id in cur_rec.children:
				if child_id not in deleted_hog_ids:
					stack.append(child_id)

	# 按层级顺序遍历树（根→叶）
	for node in tree.traverse(strategy="levelorder"):
		node_name = node.name
		logger.info('Building deletion map: processing node %s', node_name)
		node_hogs = node_to_hogs.get(node_name, [])
		# 过滤掉已标记删除的HOG
		node_hogs = [(hid, rec) for hid, rec in node_hogs if hid not in deleted_hog_ids]

		# 根节点：处理同SOG的重复HOG
		if node.is_root():
			hogs_by_og = defaultdict(list)
			for hid, rec in node_hogs:
				hogs_by_og[rec.og_id].append((hid, rec))

			for og_id, hogs in hogs_by_og.items():
				if len(hogs) > 1:
					# 按物种数降序，保留物种最多的
					hogs.sort(key=lambda x: (hog_species_count[x[0]], x[0]), reverse=True)
					for hid, rec in hogs[1:]:
						deleted_genes.update(rec.genes)
						mark_hog_and_descendants(hid)
			# 重新过滤node_hogs（因为可能有新删除的）
			node_hogs = [(hid, rec) for hid, rec in node_hogs if hid not in deleted_hog_ids]

		# 叶子节点：特殊处理，取父节点的HOGs（保持原逻辑）
		if node.is_leaf() and node.up:
			node_hogs = node_to_hogs.get(node.up.name, [])
			node_hogs = [(hid, rec) for hid, rec in node_hogs if hid not in deleted_hog_ids]

		logger.debug('Processing node %s (%d HOGs)', node_name, len(node_hogs))
		if not node_hogs:
			logger.debug('No HOGs for node %s (may have been removed by previous step)', node_name)
		
		sorted_node_hogs = sorted(node_hogs, key=lambda x: x[0])
		
		# 处理当前节点的每一个HOG
		for hid, hog_record in sorted_node_hogs:
			if hid in deleted_hog_ids:
				continue

			if node.is_leaf():
				# 叶子节点：只保留该物种的一个基因
				species_genes = [g for g in hog_record.genes
								if gene_format_o(g)[0] == node_name and g not in deleted_genes]
				if len(species_genes) > 1:
					if gene_degree:
						species_genes.sort(key=lambda g: (gene_degree.get(g, 0), g), reverse=True)
					else:
						species_genes.sort()
					deleted_genes.update(species_genes[1:])
			else:
				# 中间节点：处理子HOG，每个子节点下只保留一个子HOG
				children_by_node = defaultdict(list)
				for child_id in sorted(hog_record.children):
					child_hog = all_hogs.get(child_id)
					if child_hog and child_id not in deleted_hog_ids:
						children_by_node[child_hog.node_id].append((child_id, child_hog))

				for node_id, child_hogs in children_by_node.items():
					if len(child_hogs) > 1:
						selected_id, selected_rec = child_hogs[0]
						for child_id, child_rec in child_hogs[1:]:
							species_count = hog_species_count[child_id]
							selected_species_count = hog_species_count[selected_id]
							if species_count > selected_species_count:
								deleted_genes.update(selected_rec.genes)
								mark_hog_and_descendants(selected_id)
								selected_id = child_id
								selected_rec = child_rec
							else:
								deleted_genes.update(child_rec.genes)
								mark_hog_and_descendants(child_id)
	return deleted_genes


def process_og_with_hog(og_file, hog_args, output_file,
						restore_gene=False, restore_log=None):
	"""
	根据HOG信息处理OG文件（优化版：缓存基因物种映射，减少重复解析）

	:param og_file: 输入OG文件路径 (MCL格式)
	:param hog_args: 包含HOG参数的字典 (ogfile, orthfiles, sptreefile, paralog)
	:param output_file: 输出文件路径
	:param restore_gene: 当一个物种没有任何基因保留时，是否恢复该物种的第一个基因
	:param restore_log: 记录恢复基因的日志文件路径
	"""
	start_time = time.time()

	# 加载HOG信息（包含预计算的物种数）
	all_hogs, tree, hog_species_count = load_hog_info_from_class(
		ogfile=hog_args['ogfile'],
		orthfiles=hog_args['orthfiles'],
		sptreefile=hog_args['sptreefile'],
		paralog=hog_args.get('paralog', False)
	)
	
	load_hog_end = time.time()
	logger.info('load_hog took %.2f s', load_hog_end - start_time)

	# 构建基因degree图

	gene_degree_graph = ColinearGroups(hog_args['orthfiles'], noparalog=False).graph
	gene_degree = dict(gene_degree_graph.degree())
	
	gene_degree_end = time.time()
	logger.info('Building gene_degree took %.2f s', gene_degree_end - load_hog_end)

	# 构建待删除基因集合
	all_deleted_genes = build_deletion_map(all_hogs, tree, hog_species_count, gene_degree)
	end_time = time.time()
	logger.info('build_deletion_map took %.2f s', end_time - gene_degree_end)

	# 缓存基因->物种映射，避免重复调用 gene_format_o
	gene_to_species = {}

	# 处理OG文件
	with open(output_file, 'w') as fout, \
		(open(restore_log, 'w') if restore_log else nullcontext()) as log_file:

		line_count = 0
		for og in OrthoMCLGroup(og_file):
			og_id = og.ogid
			genes = og.genes

			# 按物种分组基因（利用缓存）
			genes_by_species = defaultdict(list)
			for gene in genes:
				sp = gene_to_species.get(gene)
				if sp is None:
					sp, _ = gene_format_o(gene)
					gene_to_species[gene] = sp
				if sp is not None:
					genes_by_species[sp].append(gene)

			kept_genes = set()
			restored_species = []

			for species, sp_genes in genes_by_species.items():
				active = [g for g in sp_genes if g not in all_deleted_genes]
				if active:
					kept_genes.update(active)
				elif sp_genes and restore_gene:
					# 恢复该物种的第一个基因
					orig = sp_genes[0]
					kept_genes.add(orig)
					restored_species.append(f"{species}:{orig}")
				elif sp_genes:
					# 未启用恢复，但仍记录日志
					restored_species.append(f"{species}:{sp_genes[0]}")

			if restored_species and log_file:
				log_file.write(f"{og_id}\t{';'.join(restored_species)}\n")

			final_genes = sorted(kept_genes) if restore_gene else sorted(kept_genes - all_deleted_genes)
			if final_genes:
				fout.write(f"{og_id}: {' '.join(final_genes)}\n")

			line_count += 1
			# 每5000行刷新一次（可选，但一般由系统自动缓冲，这里保留适度刷新）
			if line_count % 5000 == 0:
				fout.flush()
		fout.flush()


def main():
	parser = argparse.ArgumentParser(description="从HOG信息创建单拷贝OG文件")
	parser.add_argument("-og", "--ogfile", required=True, help="输入的OG文件（MCL格式）")
	parser.add_argument("-hog_og", "--hog_ogfile", required=True, help="HOG所需的OG文件（与-og相同或不同）")
	parser.add_argument("-orth", "--orthfiles", nargs='+', required=True, help="正交文件列表（空格分隔）")
	parser.add_argument("-tree", "--sptreefile", required=True, help="物种树Newick文件")
	parser.add_argument("-o", "--output", required=True, help="输出文件路径")
	parser.add_argument("--paralog", action="store_true", help="是否包含旁系同源")
	parser.add_argument("--restore_gene", action="store_true", help="当一个物种没有任何基因保留时，恢复该物种的第一个基因")
	parser.add_argument("--restore_log", help="记录恢复基因的日志文件")

	args = parser.parse_args()

	hog_args = {
		'ogfile': args.hog_ogfile,
		'orthfiles': args.orthfiles,
		'sptreefile': args.sptreefile,
		'paralog': args.paralog
	}

	process_og_with_hog(
		og_file=args.ogfile,
		hog_args=hog_args,
		output_file=args.output,
		restore_gene=args.restore_gene,
		restore_log=args.restore_log,
	)


if __name__ == "__main__":
	main()