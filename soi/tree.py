import re
from ete3 import Tree

class SpeciesTree:
	pass

def number_nodes(sptreefile):
	treestr = convertNHX(sptreefile)
#	print(treestr)
	tree = load_tree_smart(treestr)
#	try: tree = Tree(treestr, format=0)
#	except :
#		tree = Tree(treestr, format=1)
	i = 0
	for node in tree.traverse():
		node.show = True
		if node.is_leaf():
			continue
		if node.name:
			continue
		name = 'N{}'.format(i)
		node.name = name
		i += 1
	# 压缩单子节点：删掉只有一个孩子的内部节点，祖孙直连
	for node in list(tree.traverse()):
		if not node.is_leaf() and not node.is_root() and len(node.children) == 1:
			node.delete(prevent_nondicotomic=False)
	tree.write(outfile=sptreefile + ".labeled.nwk", format=1)
	return tree 


def load_tree_smart(treestr):
    """
    智能判断：检测内部节点标签是纯数字(support)还是文本(name)
    """
    # 提取所有 ")LABEL:" 或 ")LABEL;" 模式的内部节点标签
    # 匹配 ) 后到 : 或 ; 或 ) 前的内容
    internal_labels = re.findall(r'\)([^:;,\(\)]+)(?=:)', treestr)
    
    # 判断是否有非数字标签（允许空，允许小数）
    has_text_label = any(
        label and not re.match(r'^\d*\.?\d+$', label) 
        for label in internal_labels
    )
    
    if has_text_label:
        # 有文本标签（如 N1），用 format=1
        return Tree(treestr, format=1)
    else:
        # 只有数字或空，用 format=0（保留 support 值）
        return Tree(treestr, format=0)

def convert_newick(line: str) -> str:
    '''(A,B[p=2]); (A,B[p=2]:0.1); (A,B:0.1[p=2]); (A,B[&&NHX:p=2])'''
    # 1. [tag]:length → :length[&&NHX:tag]（带安全检查）
    def patch_tag_length(m):
        tag, length = m.groups()
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f":{length}[{tag}]"
    
    line = re.sub(r'\[([^]]+)\]:([\d.eE-]+)', patch_tag_length, line)

    # 2. :length[tag] → :length[&&NHX:tag]
    def patch_length_tag(m):
        length, tag = m.groups()
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f":{length}[{tag}]"
    
    line = re.sub(r':([\d.eE-]+)\[([^]]+)\]', patch_length_tag, line)

    # 3. [tag]（孤立无长度）→ [&&NHX:tag]
    def patch_isolated(m):
        tag = m.group(1)
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f"[{tag}]"
    
    line = re.sub(r'\[([^]]+)\](?![\d:\[])', patch_isolated, line)
#    print(line)
    return line
	
def convertNHX(inNwk, ):
	nwk = []
	for line in open(inNwk):
		nwk += [convert_newick(line)]
	return ''.join(nwk)


class TreePlotter:
	"""Matplotlib orthogonal tree plotter, Dsuite-Fbranch style.

	Every node (internal + leaves) sits on its own row; row order comes from
	`row_names` (Dsuite Fbranch mode: preorder over ALL nodes when
	`internal_node_to_present=True`). Internal nodes get dotted alignment
	lines to the right edge, like Fbranch's guides between tree and matrix.
	Draws into a caller-provided matplotlib Axes.

	Usage:
		tp = TreePlotter(tree, row_names, internal_node_to_present=True)
		tp.plot(ax)
	"""
	def __init__(self, tree, row_names=None, internal_node_to_present=False,
				 use_distances=False, color='k', lw=1.0,
				 dotted_lw=0.8, dotted_color='0.6', xpad=0.2,
				 leaf_dotted_color=None):
		"""tree: ete3 Tree (nodes named, e.g. via number_nodes()).
		row_names: row order list; None = auto (preorder all nodes when
		internal_node_to_present else leaf names).
		use_distances: x = time from root (branch lengths) instead of depth.
		"""
		self.tree = tree
		self.internal_node_to_present = internal_node_to_present
		self.use_distances = use_distances
		self.color = color
		self.lw = lw
		self.dotted_lw = dotted_lw
		self.dotted_color = dotted_color
		self.leaf_dotted_color = leaf_dotted_color or dotted_color
		self.xpad = xpad
		if row_names is None:
			if internal_node_to_present:
				row_names = [n.name for n in tree.traverse('preorder')]
			else:
				row_names = tree.get_leaf_names()
		self.row_names = list(row_names)
		self._layout()

	def _layout(self):
		"""Assign x (depth/time) and y (row index) to every node."""
		tree = self.tree
		row_idx = {}
		for i, nm in enumerate(self.row_names):
			if nm not in row_idx:
				row_idx[nm] = i
		max_dist = tree.get_farthest_leaf()[1] if self.use_distances else 0.0
		ypos = {}
		xpos = {}
		def layout(node, depth):
			xpos[node] = (max_dist - tree.get_distance(node)) if self.use_distances else depth
			# recurse first (children must be placed before the fallback mean)
			for ch in node.children:
				layout(ch, depth + 1)
			if node.name in row_idx:
				ypos[node] = row_idx[node.name]
			else:
				# name not in row list (e.g. leaves-only mode with unlisted
				# leaves): fall back to mean of children rows
				ys = [ypos[ch] for ch in node.children if ch in ypos]
				ypos[node] = sum(ys) / len(ys) if ys else None
		layout(tree, 0)
		self.ypos = ypos
		self.xpos = xpos

	def plot(self, ax, show_dotted=True):
		"""Draw orthogonal tree into ax; dotted guides from internal nodes."""
		tree = self.tree
		ypos, xpos = self.ypos, self.xpos
		for node in tree.traverse():
			if node.is_root():
				continue
			yp = ypos.get(node)
			if yp is None:
				continue
			par = node.up
			x, xp = xpos[node], xpos[par]
			# horizontal branch at child y
			ax.plot([xp, x], [yp, yp], '-', color=self.color, lw=self.lw,
					zorder=1, solid_capstyle='butt')
			# vertical connector at parent x
			y_par = ypos.get(par)
			if y_par is not None:
				ax.plot([xp, xp], [y_par, yp], '-', color=self.color,
						lw=self.lw, zorder=1, solid_capstyle='butt')
		# root stub
		y_root = ypos.get(tree)
		x_root = xpos[tree]
		if y_root is not None:
			ax.plot([x_root - self.xpad, x_root], [y_root, y_root], '-',
					color=self.color, lw=self.lw, zorder=1)
		# dotted alignment guides from internal nodes to the right edge
		xmax = max(xpos.values())
		if show_dotted and self.internal_node_to_present:
			for node in tree.traverse():
				if node.is_root():
					continue
				yp = ypos.get(node)
				if yp is None:
					continue
				# leaf guides start at the branch tip; internal guides start at the node
				color = self.leaf_dotted_color if node.is_leaf() else self.dotted_color
				ax.plot([xpos[node], xmax + self.xpad], [yp, yp], ':',
						color=color, lw=self.dotted_lw, zorder=0)
		# axes: rows span the y range, top-down when row mode
		if self.internal_node_to_present:
			ax.set_ylim(-0.5, len(self.row_names) - 0.5)
			ax.invert_yaxis()
		else:
			ys = [v for v in ypos.values() if v is not None]
			ax.set_ylim(min(ys) - 0.5, max(ys) + 0.5)
		ax.set_xlim(-self.xpad, xmax + (self.xpad if show_dotted and self.internal_node_to_present else 0))
		ax.axis('off')