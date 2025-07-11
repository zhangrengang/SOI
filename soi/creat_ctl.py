import sys
import re
#from mcscan import Collinearity
from .small_tools import parse_kargs
def main(outCtl = sys.stdout, pngsize=1500):
	subcmd = sys.argv[1]
	if subcmd == 'by_match':
		collinearity, refSP, qrySP = sys.argv[2:5]
		chrAlist, chrBlist = collinearity_counter(collinearity, refSP, qrySP)
	elif subcmd == 'by_genes':
		inChrList, refSP, qrySP = sys.argv[2:5]
		chrAlist, chrBlist = filter_chrs(inChrList, refSP, qrySP)
	elif subcmd == 'jcvi':
		kargs = parse_kargs(sys.argv)
		collinearity, cfg = sys.argv[2:4]
		jcvi_seqids(collinearity, cfg, outCtl=outCtl, **kargs)
		return
	else:
		raise ValueError('Unknown sub command: {}'.format(subcmd))
	print('%s\n%s' % (pngsize, pngsize), file=outCtl)
	print(','.join(chrBlist), file=outCtl)
	print(','.join(chrAlist), file=outCtl)
def jcvi_seqids(collinearity, cfg, ref=None, outCtl=sys.stdout, by='match'):
	sps = [line.strip().split()[0] for line in open(cfg)]
	if not ref:
		ref = sps[0]
	chrList = []
	here = False
	for i, sp in enumerate(sps):
		if sp == ref:
			try: chrList += [chrAlist]
			except: here = True
			continue
		chrAlist, chrBlist = collinearity_counter(collinearity, ref, sp) if by == 'match' else \
			filter_chrs(collinearity, ref, sp)
		if here:
			chrList += [chrAlist]
			here = False
		chrList += [chrBlist]
	for chrBlist in chrList:
		print(','.join(chrBlist), file=outCtl)

def filter_chrs(inChrList, refSP, qrySP, min_genes=200, 
		ref_min_genes=None, qry_min_genes=None,
		min_ratio=0.5):
	if ref_min_genes is None:
		ref_min_genes = min_genes
	if qry_min_genes is None:
		qry_min_genes = min_genes
	chrAlist, chrBlist = [],[]
	ref_total_genes, ref_filtered_genes = 0, 0
	qry_total_genes, qry_filtered_genes = 0, 0
	for line in open(inChrList):
		temp = line.strip().split('\t')
		chr, CHR, SP, geneN = temp[:4]
		geneN = int(geneN)
#		total_genes += geneN
#		if geneN < min_genes:
#			continue
#		filtered_genes += geneN
#		if SP == qrySP:
#			print >>sys.stderr, SP,qrySP,chr, CHR, is_chr0(CHR), is_chr0(chr, raw=False)
		if is_chr0(CHR) or is_chr0(chr, raw=False):
			continue
		pattern = '\w{1,2}\d+'
		if not re.compile(pattern).match(chr):
			continue
		pattern = '{}\d+'.format(refSP)
		if SP == refSP or re.compile(pattern).match(chr):
			ref_total_genes += geneN
			if geneN >= ref_min_genes:
				chrAlist += [chr]
				ref_filtered_genes += geneN
		#if not (SP == refSP or SP == qrySP):
		#	continue
		pattern = '{}\d+'.format(qrySP)
#		print >>sys.stderr, SP,qrySP,chr
		if SP == qrySP or re.compile(pattern).match(chr):
	#		print >>sys.stderr, temp
			qry_total_genes += geneN
			if geneN >= qry_min_genes:
				chrBlist += [chr]
				qry_filtered_genes += geneN
#	print >>sys.stderr, chrAlist, chrBlist, refSP, qrySP
	ref_ratio = 1.0*ref_filtered_genes/ref_total_genes
	qry_ratio = 1.0*qry_filtered_genes/qry_total_genes
	chrAlist = sort_version(chrAlist)
	chrBlist = sort_version(chrBlist)
	#if not (chrAlist and chrBlist): # too fragment
	min_genes = 10
	if ref_ratio < min_ratio and qry_ratio<min_ratio:
		return filter_chrs(inChrList, refSP, qrySP, min_genes=min_genes)
	elif ref_ratio < min_ratio:
		return filter_chrs(inChrList, refSP, qrySP, ref_min_genes=min_genes)
	elif qry_ratio < min_ratio:
		return filter_chrs(inChrList, refSP, qrySP, qry_min_genes=min_genes)
	return chrAlist, chrBlist
def get_good_chrs(inChrList, min_genes=200, min_ratio=0.8):
	good_chrs = []
	N, n  = 0, 0
	for line in open(inChrList):
		temp = line.strip().split('\t')
		chr, CHR, SP, geneN = temp[:4]
		geneN = int(geneN)
		if is_chr0(CHR) or is_chr0(chr, raw=False):
			continue
		N += geneN
		pattern = '\w{1,2}\d+'
		if not re.compile(pattern).match(chr):
			continue
		if geneN < min_genes:
			continue
		n += geneN
		good_chrs += [chr]
	if N == 0:
		raise ValueError('No genes detected in `{}`. Please check it'.format(inChrList))
	ratio = 1.0* n /N
	if ratio < min_ratio:
		return get_good_chrs(inChrList, min_genes=min_genes-25,min_ratio=min_ratio)
	return good_chrs
	
class Chrs:
	def __init__(self, inChrList):
		self.inChrList = inChrList
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.inChrList):
			yield ChrLine(line)
class ChrLine:
	def __init__(self, line):
		temp = line.strip().split('\t')
		chr, CHR, SP, geneN = temp[:4]
		geneN = int(geneN)
		self.chr, self.CHR, self.SP, self.geneN = chr, CHR, SP, geneN
def is_chr0(chr, raw=True):	# not a ture chromosome
#	print >>sys.stderr, chr
	match = re.compile(r'(\d+)').search(chr)
#	if match:
#		print >>sys.stderr, chr, match.groups()
	if match and match.groups()[0] == '0':	# chr0
		return True
	elif raw and re.compile(r'[MP]t', re.I).match(chr) and not re.compile(r'ptg', re.I).match(chr) :	# Mt
		return True
	elif not re.compile(r'\d', re.I).search(chr):	# chrUn
		return True
	else:
		return False
def collinearity_counter(collinearity, refSP, qrySP, sort_ref=True, sort_qry=False, min_genes=100):
	from .mcscan import Collinearity
	d_chrB = {}
	d_count = {}
	for rc in Collinearity(collinearity):
		if not rc.species == (refSP, qrySP): # ( rc.species == (refSP, qrySP) or rc.species == (qrySP, refSP) ):
			continue
		chrA, chrB = rc.chrs 
#		print chrA, chrB, rc.species
		if (rc.species1, rc.species2 ) == (qrySP, refSP):
			chrA, chrB = chrB, chrA # A: ref, B:qry
		try: d_chrB[chrB] += [(rc.score, chrA)]
		except KeyError: d_chrB[chrB] = [(rc.score, chrA)]
		for chr in rc.chrs:
			try: d_count[chr] += rc.N
			except KeyError: d_count[chr] = rc.N

#	print d_chrB	
#	print d_count
	# best match
	d_chrAB = {}
	for chrB, info in list(d_chrB.items()):
		score, chrA = sorted(info)[-1]
		try: d_chrAB[chrA].append(chrB)
		except KeyError: d_chrAB[chrA] = [chrB]

	chrAlist, chrBlist = [],[]
	for chrA in sort_version(list(d_chrAB.keys())):
		chrBs = d_chrAB[chrA]
		chrAlist.append(chrA)
		chrBlist += chrBs
	if sort_qry:
		chrBlist = sort_version(chrBlist)
	
	chrAlist = [chr for chr in chrAlist if d_count[chr] >= min_genes]
	chrBlist = [chr for chr in chrBlist if d_count[chr] >= min_genes]
	return chrAlist, chrBlist
def sort_version(alist):
	blist = []
	for v in alist:
		try:
			i = int(re.compile(r'(\d+)').search(v).groups()[0])
		except AttributeError:
			i = 0
		try:
			p = re.compile(r'([^\d]*)').match(v).groups()[0]
		except AttributeError:
			p = ''
		blist.append((p, i))
	return [y for x,y in sorted(zip(blist,alist))]

if __name__ == '__main__':
	main()

