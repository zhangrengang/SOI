import sys
import os
import argparse
from .RunCmdsMP import logger
from .__version__ import version

bindir = os.path.dirname(os.path.realpath(__file__))


def args_common(parser):
	pass


def args_dotplot(parser):
	from .dot_plotter import dotplot_args
	args = dotplot_args(parser)
	return dict(args=args)


def func_dotplot(**kargs):
	from .dot_plotter import xmain as dot_plotter
	dot_plotter(**kargs)


def args_filter(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files output from MCscanX, WGDI, or MCscan/JCVI. \
[required]")
	parser.add_argument('-o', '-orthology', required=True,  type=str,  nargs='*',
						dest='orthologs',  metavar='FOLDER/FILE',
						help="Orthologues output from OrthoFinder (folder), or OrthoMCL (file). \
[required]")
	parser.add_argument('-c', '-cutoff',  type=float, default=0.6,
						dest='min_ratio',  metavar='FLOAT',
						help="Cutoff (lower limit) of Orthology Index (rataining blocks > this cutoff) [default=%(default)s]")
	parser.add_argument('-u', '-upper',  type=float, default=1,
						dest='max_ratio',  metavar='FLOAT',
						help="Upper limit of Orthology Index (rataining blocks <= this cutoff) [default=%(default)s]")
	parser.add_argument('-n', '-min_n',  type=int, default=0,
						dest='min_n',  metavar='INT',
						help="Minimum gene number in a block [default=%(default)s]")
	parser.add_argument('-g', '-gff', type=str,
						dest='gff',  metavar='FILE',
						help="Gff file. [required for `-d`]")
	parser.add_argument('-d', '-min_dist',  type=int, default=None,
						dest='min_dist',  metavar='INT',
						help="Minimum distance to remove a tandem repeated block [default=%(default)s]")
	parser.add_argument('-stat', default=None,
						dest='out_stats', type=str,
						help="Output stats by species pairs. [default=%(default)s]")
	parser.add_argument('-oo', default=False,
						dest='output_orthology', action='store_true',
						help="Output retained orthology instead of synteny. [default=%(default)s]")


def func_filter(**kargs):
	from .mcscan import identify_orthologous_blocks
	identify_orthologous_blocks(**kargs)


def args_cluster(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files from `filter` sub-command. [required]")
	parser.add_argument('-o', '-orthology', type=str,  nargs='*',
						dest='orthologs',  metavar='FOLDER/FILE', default=None,
						help="Orthologues output from OrthoFinder (folder), or OrthoMCL (file). \
This will use Orthology Index as weight for MCL [default=%(default)s]")
	parser.add_argument('-I', '-inflation', type=float, default=1.5, metavar='FLOAT',
						dest='inflation',
						help="Inflation for MCL (varying this parameter affects granularity) \
[default=%(default)s]")
	parser.add_argument('-outgroup', type=str, default=None, metavar='TAXON/FILE',
						dest='outgroup', nargs='*',
						help="Outgroups to exclude from orthogroups (prior to `-ingroup`) \
[default=%(default)s]")
	parser.add_argument('-ingroup', type=str, default=None, metavar='TAXON/FILE',
						dest='ingroup', nargs='*',
						help="Ingroups that are only included [default=%(default)s]")
	parser.add_argument('-prefix', type=str, default='cluster',
						dest='outpre',
						help="Output prefix [default=%(default)s]")


def func_cluster(**kargs):
	from .mcscan import cluster_by_mcl
	cluster_by_mcl(**kargs)


def args_outgroup(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files from `filter` sub-command. [required]")
	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='orthogroup',  metavar='FILE',
						help="Orthogroups output from `cluster` sub-command. [required]")
	parser.add_argument('-outgroup', type=str, required=True, metavar='TAXON',
						dest='outgroup', nargs='*',
						help="Outgroups to include to orthogroups [required]")
	parser.add_argument('-cutoff',  type=float, default=0.2,
						dest='min_ratio',  metavar='FLOAT',
						help="Cutoff (lower limit) of links to outgroup genes [default=%(default)s]")


def func_outgroup(**kargs):
	from .mcscan import cluster_add_outgroup
	cluster_add_outgroup(**kargs)


def args_phylo_common(parser):
	parser.add_argument('-mc', '-max_copies', type=float, default=6,
						dest='max_copies', metavar='INT',
						help="To limit a common maximum copy number for every species. [default=%(default)s]")
	parser.add_argument('-sc', '-singlecopy', default=None,
						dest='singlecopy', action='store_true',
						help="Only retrieve singlecopy genes (=`-max_copies 1`). [default=%(default)s]")
	parser.add_argument('-ss', '-spsd', type=str, default=None,
						dest='spsd',  metavar='FILE',
						help="To limit a specific copy number for each species (format: 'TAXON<tab>NUMBER'). [default=%(default)s]")
	parser.add_argument('-fmt', type=str, default='orthomcl',
						dest='source', choices=['orthomcl', 'orthofinder', 'mcscanx'],
						help="Format of `-orthogroup` input. [default=%(default)s]")


def args_stats(parser):
	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='input',  metavar='FILE',
						help="Orthogroups output from `cluster` or `outgroup` sub-commands. [required]")
	parser.add_argument('-mm', '-max_missing', type=float, default=0.4,
						dest='max_taxa_missing', metavar='FLOAT',
						help="To allow maximum ratio of missing species. [default=%(default)s]")
	args_phylo_common(parser)


def func_stats(**kargs):
	from .mcscan import orthomcl_stats
	orthomcl_stats(**kargs)


def args_phylo(parser):
	import uuid
	uid = uuid.uuid1()
	default_tmpdir = './tmp-{}'.format(uid)

	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='input',  metavar='FILE',
						help="Orthogroups output from `cluster` or `outgroup` sub-commands. \
[required]")
	parser.add_argument('-pep', required=True,  type=str,
						dest='pep',  metavar='PEP FILE',
						help="Protein fasta file. [required]")
	parser.add_argument('-cds', type=str, default=None,
						dest='cds',  metavar='CDS FILE',
						help="CDS fasta file. [default=%(default)s]")
	parser.add_argument('-both', default=False,
						dest='both', action='store_true',
						help="To use both CDS and PEP to build gene trees (only valid when `-cds` is true). \
[default: %(default)s]")
	parser.add_argument('-root', '-outgroup', type=str, metavar='TAXON',
						dest='root', nargs='*', default=None,
						help="Outgroups to root gene trees [default=%(default)s]")
	parser.add_argument('-pre', '-prefix', type=str, default='sog',
						dest='suffix', metavar='STR',
						help="Output prefix. [default=%(default)s]")
	parser.add_argument('-mm', '-max_missing', type=float, default=0.4,
						dest='max_taxa_missing', metavar='FLOAT',
						help="To allow maximum ratio of missing species. [default=%(default)s]")
	args_phylo_common(parser)

	parser.add_argument('-aligner', type=str, default='muscle',
						dest='aligner', metavar='STR', choices=['muscle', 'mafft'], 
						help="Aligner: muscle (v5/v3) or mafft. [default=%(default)s]")

	parser.add_argument('-only_aln', default=False,
						dest='onlyaln', action='store_true',
						help="Only aligning sequences, to skip trimal and iqtree. [default: %(default)s]")
	parser.add_argument('-concat', default=False,
						dest='concat', action='store_true',
						help="To concatenate alignments for tools such as IQTREE \
(valid when `-singlecopy` is true). [default: %(default)s]")
	parser.add_argument('-trimal_opts', type=str, default='-automated1',
						dest='trimal_opts',  metavar='STR',
						help="TrimAl options. [default='%(default)s']")
	parser.add_argument('-iqtree_opts', type=str, default='',
						dest='iqtree_opts',  metavar='STR',
						help="IQ-TREE options. [default='%(default)s']")
	parser.add_argument('-fast', default=False,
						dest='fast', action='store_true',
						help="Speedup IQ-TREE by restricting model set (JTT for PEP and GTR for CDS). \
[default: %(default)s]")
	parser.add_argument('-p', '-ncpu', type=int, default=20,
						dest='ncpu', metavar='INT',
						help="Number of processors. [default=%(default)s]")
	parser.add_argument('-tmp', '-tmpdir', type=str, default=default_tmpdir,
						dest='tmpdir',  metavar='FOLDER',
						help="Temporary folder. [default=%(default)s]")
	parser.add_argument('-clean', default=False,
						dest='clean', action='store_true',
						help="Cleanup temporary folder. [default: %(default)s]")
	# overwrite


def func_phylo(**kargs):
	from .mcscan import orthomcl_to_astral
	orthomcl_to_astral(**kargs)


def makeArgs():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Play with Orthology Index',
	)
	parser.add_argument(
		'-v', '--version',
		action='version',
		version=version
	)
	# subcommands
	subparsers = parser.add_subparsers(help='sub-command help')
	parser_dot = subparsers.add_parser('dotplot',
									   help='Generate Ks/OI-colored dot plots')
	args_dotplot(parser_dot)
	parser_flt = subparsers.add_parser('filter',
									   help='Filter synteny with Orthology Index (standard output)')
	args_filter(parser_flt)
	parser_clst = subparsers.add_parser('cluster',
										help='Cluster syntenic orthogroups (SOGs)')
	args_cluster(parser_clst)
	parser_ogrp = subparsers.add_parser('outgroup',
										help='Add outgroups for SOGs from synteny')
	args_outgroup(parser_ogrp)
	parser_phylo = subparsers.add_parser('phylo',
										 help='Build gene trees from SOGs')
	args_phylo(parser_phylo)
	parser_stats = subparsers.add_parser('stats',
										 help='Make statistics of SOGs for phylogeny')
	args_stats(parser_stats)

	if len(sys.argv) == 1:
		parser.print_help(sys.stderr) 
		sys.exit(1) 
	args = parser.parse_args()
	return args


FUNC = {
	'dotplot': func_dotplot,
	'filter': func_filter,
	'cluster': func_cluster,
	'outgroup': func_outgroup,
	'phylo': func_phylo,
	'stats': func_stats,
}


def main():
	args = makeArgs()  # options
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	key = sys.argv[1]
	func = FUNC[key]  # functions
	func(**args.__dict__)  # execute
	logger.info('Completed\n')


if __name__ == '__main__':
	main()
