#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
takr_colored_graph.py - ColoredGraph: 彩色邻接图 + 环检测 + 路径覆盖

每条边记录一组颜色标签 (child_id, chromosome_id)，
一条边可以有多个颜色（多个孩子共享该邻接）。
颜色集为空时自动移除边。

使用无向图 (nx.Graph) 存储，方向在路径覆盖后由 AncestralAdjacencyGraph 处理。
"""

import logging
import math
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

from .takr_events import TAKREvent
from . import style  # global mpl settings (no-op if already imported)

logger = logging.getLogger(__name__)


def _merge_chromosome_paths(path_a, path_b):
    """Merge two subgenome chromosome paths into one pre-WGD path.
    
    Strategy: union of HOGs, keep path_a order first, append unique from path_b.
    """
    seen = set()
    merged = []
    for h in path_a:
        if h not in seen:
            merged.append(h)
            seen.add(h)
    for h in path_b:
        if h not in seen:
            merged.append(h)
            seen.add(h)
    return merged


def _assign_pair_colors(labels):
    """用 golden ratio hash 为 label 列表分配 HSL 色板。

    Returns: {label: '#RRGGBB'}
    """
    if not labels:
        return {}
    result = {}
    for i, label in enumerate(labels):
        hue = (i * 0.618033988749895) % 1.0
        r, g, b = _hsl_to_rgb(hue, 0.55, 0.50)
        result[label] = f'#{r:02x}{g:02x}{b:02x}'
    return result


def _hsl_to_rgb(h, s, l):
    """HSL (h,s,l ∈ [0,1]) → (r,g,b ∈ [0,255])."""
    if s == 0:
        v = int(l * 255)
        return (v, v, v)

    def _hue2rgb(p, q, t):
        if t < 0:
            t += 1
        if t > 1:
            t -= 1
        if t < 1 / 6:
            return p + (q - p) * 6 * t
        if t < 1 / 2:
            return q
        if t < 2 / 3:
            return p + (q - p) * (2 / 3 - t) * 6
        return p

    q = l * (1 + s) if l < 0.5 else l + s - l * s
    p = 2 * l - q
    r = int(_hue2rgb(p, q, h + 1 / 3) * 255)
    g = int(_hue2rgb(p, q, h) * 255)
    b = int(_hue2rgb(p, q, h - 1 / 3) * 255)
    return (r, g, b)


class ColoredGraph(nx.DiGraph):
    """彩色邻接图 — 一条边可有多色 (child_id, chrom_idx)。"""

    def __init__(self, hog_level: str = ''):
        super().__init__()  # nx.DiGraph
        self.hog_level = hog_level
        self.events: List[TAKREvent] = []
        # 记录每个 child 的端粒集，用于 indel 检测和路径覆盖
        self._child_telomeres: Dict[str, Set] = {}
        # 记录每个 child 的所有 HOG 集合
        self._child_hogs: Dict[str, set] = {}
        # 记录每个 child 的染色体数，用于验证
        self._child_chrom_counts: Dict[str, int] = {}
        # 记录每个 child 的染色体路径，用于单孩子块压缩
        self._child_chromosomes: Dict[str, list] = {}

    def _neighbors(self, node):
        """有向图邻居 = 前驱 ∪ 后继。"""
        return set(self.predecessors(node)) | set(self.successors(node))

    def _degree(self, node):
        """有向图度数 = in_degree + out_degree。"""
        return self.in_degree(node) + self.out_degree(node)

    # ==================== 构建 ====================

    def _block_edge_data(self, b1, b2):
        """获取块间边数据（有向图：自动找 b1→b2 或 b2→b1）。"""
        bg = self._block_graph
        if bg.has_edge(b1, b2):
            return bg[b1][b2]
        if bg.has_edge(b2, b1):
            return bg[b2][b1]
        return {}

    def _block_has_edge(self, b1, b2) -> bool:
        """块间边是否存在（任意方向）。"""
        bg = self._block_graph
        return bg.has_edge(b1, b2) or bg.has_edge(b2, b1)

    def _block_add_edge(self, b1, b2, **attrs):
        """添加块间有向边 b1→b2（保留方向）。"""
        self._block_graph.add_edge(b1, b2, **attrs)

    def _block_remove_edge(self, b1, b2):
        """移除块间边（任意方向）。"""
        bg = self._block_graph
        if bg.has_edge(b1, b2):
            bg.remove_edge(b1, b2)
        elif bg.has_edge(b2, b1):
            bg.remove_edge(b2, b1)

    def add_synteny_edge(self, h1, h2, child_id: str, chrom_idx: int, direction: int = 1):
        """添加有向边 h1→h2（direction=+1）或 h2→h1（direction=-1）。

        有向图: 边方向由 direction 决定。
        """
        if direction == -1:
            h1, h2 = h2, h1
        key = (h1, h2)

        if self.has_edge(key[0], key[1]):
            colors = self[key[0]][key[1]]['colors']
            colors.add((child_id, chrom_idx))
            directions = self[key[0]][key[1]].get('directions', set())
            directions.add((child_id, chrom_idx, 1))
            self[key[0]][key[1]]['directions'] = directions
        else:
            self.add_edge(key[0], key[1],
                                 colors={(child_id, chrom_idx)},
                                 directions={(child_id, chrom_idx, 1)})

    def add_child(self, source_id: str, child_graph) -> int:
        """把一个孩子的所有染色体邻接以该孩子的颜色加入。

        方向: 孩子中 HOG 的顺序决定方向。
        无向边 {A, B} 的方向是相对于规范顺序 (A, B where A < B):
          - 孩子有 A→B: direction = +1
          - 孩子有 B→A: direction = -1
        """
        chrom_count = 0
        child_hogs = set()
        child_tels = set()
        child_chrom_paths = []

        # 优先用 chrom_hogs（保留重复节点），回退到 chromosomes
        chrom_source = getattr(child_graph, 'chrom_hogs', None)
        if chrom_source:
            chrom_iter = sorted(chrom_source.items())
        else:
            chrom_iter = enumerate(child_graph.chromosomes)

        for chrom_idx, chrom_nodes in chrom_iter:
            chrom_count += 1
            hogs = [n for n in chrom_nodes if n not in child_graph.telomeres]
            child_hogs.update(hogs)
            child_chrom_paths.append(hogs)
            for i, h in enumerate(hogs):
                # 节点属性: sources, telomere（per-child 集合，GFA/算法共用）
                self.add_node(h)  # 确保节点存在
                self.nodes[h].setdefault('sources', set()).add(
                    (source_id, chrom_idx))
                if i < len(hogs) - 1:
                    self.add_synteny_edge(h, hogs[i + 1], source_id, chrom_idx,
                                  direction=1)

        # 端粒 = 每条染色体的首尾 HOG
        for hogs in child_chrom_paths:
            if hogs:
                for h in (hogs[0], hogs[-1]):
                    if self.has_node(h):
                        self.nodes[h].setdefault('telomere', set()).add(source_id)
                        child_tels.add(h)

        self._child_telomeres[source_id] = child_tels
        self._child_hogs[source_id] = child_hogs
        self._child_chrom_counts[source_id] = chrom_count
        self._child_chromosomes[source_id] = child_chrom_paths
        return chrom_count

    # ==================== 查询 ====================

    def get_directions(self, h1, h2) -> set:
        """获取边 (h1, h2) 的方向集合: {(child_id, chrom_idx, direction), ...}"""
        if not self.has_edge(h1, h2):
            return set()
        return self[h1][h2].get('directions', set())

    def edge_has_direction_conflict(self, h1, h2) -> bool:
        """检查边 (h1, h2) 是否有方向冲突 — 不同孩子方向相反。

        如果 child A 认为 h1→h2 (+1)，child B 认为 h2→h1 (-1)，
        则存在方向冲突，说明发生了 inversion。
        """
        directions = self.get_directions(h1, h2)
        if len(directions) <= 1:
            return False
        has_forward = any(d == 1 for _, _, d in directions)
        has_reverse = any(d == -1 for _, _, d in directions)
        return has_forward and has_reverse

    def find_indel_shortcuts(self, max_span: int = 20) -> List[Tuple]:
        """找 indel shortcuts: 跨越边 (spanning edge) 且无方向冲突。

        关键区分 indel vs 重排:
        - indel: 中间 HOGs 在另一个孩子中不存在（真正的插入/删除）
        - 重排: 中间 HOGs 在另一个孩子中存在但顺序不同（inversion/RT）

        Returns: [(h1, h2, child_id, spanned_hogs), ...]
        """
        shortcuts = []
        n_unique = 0
        n_dir_conflict = 0
        n_not_in_other = 0
        n_no_path = 0
        n_short_path = 0
        n_too_long = 0
        n_rearrangement = 0
        for h1, h2 in self.unique_edges():
            n_unique += 1
            if self.edge_has_direction_conflict(h1, h2):
                n_dir_conflict += 1
                continue
            colors = self.get_colors(h1, h2)
            if not colors:
                continue
            child_id, _ = next(iter(colors))
            for other_id in self.children():
                if other_id == child_id:
                    continue
                other_hogs = self._child_hogs.get(other_id, set())
                if h1 not in other_hogs or h2 not in other_hogs:
                    n_not_in_other += 1
                    continue
                path = self._shortest_path_through_hogs(h1, h2, other_id)
                if not path:
                    n_no_path += 1
                    continue
                if len(path) <= 2:
                    n_short_path += 1
                    continue
                spanned = path[1:-1]
                if len(spanned) > max_span:
                    n_too_long += 1
                    continue
                # 关键检查: 中间 HOGs 是否在 child_id 中也存在？
                # 如果存在 → 重排（不是 indel）
                child_hogs = self._child_hogs.get(child_id, set())
                n_in_child = sum(1 for h in spanned if h in child_hogs)
                if n_in_child > len(spanned) * 0.5:
                    # 超过一半的中间 HOGs 在 child_id 中也存在 → 重排
                    n_rearrangement += 1
                    continue
                shortcuts.append((h1, h2, child_id, spanned))
                break
        logger.debug("  [indel] find_indel_shortcuts: %d unique edges, "
                     "dir_conflict=%d, not_in_other=%d, no_path=%d, "
                     "short_path=%d, too_long=%d, rearrangement=%d, shortcuts=%d",
                     n_unique, n_dir_conflict, n_not_in_other, n_no_path,
                     n_short_path, n_too_long, n_rearrangement, len(shortcuts))
        return shortcuts

    def consensus_telomeres(self, min_children: int = 2) -> Set:
        """共识端粒：在 ≥min_children 个孩子中都是端粒邻接的 HOG。

        生物学原理：端粒位置保守，不会无缘无故获得或丢失。
        如果一个 HOG 在多个孩子中都与端粒相邻，它在祖先中也是端粒邻接的。
        """
        return {n for n, d in self.nodes(data=True)
                if len(d.get('telomere', set())) >= min_children}

    def child_telomere_set(self) -> Set:
        """所有孩子的端粒邻接 HOG 的并集。"""
        return {n for n, d in self.nodes(data=True) if d.get('telomere')}

    def is_telomere_adjacent(self, hog) -> bool:
        """该 HOG 是否在任意孩子中与端粒相邻。"""
        return bool(self.nodes[hog].get('telomere')) if self.has_node(hog) else False

    def telomere_preserving_color(self, cycle_edges, edge_colors) -> Optional[str]:
        """在交替色环中，找出保留端粒邻接的颜色。

        对每种颜色，检查该颜色的边是否涉及端粒 HOG。
        保留端粒的颜色 = 祖先态（不应被移除）。
        """
        # 统计每种颜色涉及的端粒 HOG 数
        color_tel_count = defaultdict(int)
        color_total = defaultdict(int)
        for (h1, h2), colors in zip(cycle_edges, edge_colors):
            for cid, _ in colors:
                color_total[cid] += 1
                if self.is_telomere_adjacent(h1) or self.is_telomere_adjacent(h2):
                    color_tel_count[cid] += 1

        # 找端粒比例最高的颜色
        best_color = None
        best_ratio = -1
        for cid in color_total:
            ratio = color_tel_count.get(cid, 0) / color_total[cid]
            if ratio > best_ratio:
                best_ratio = ratio
                best_color = cid

        # 只有当端粒比例有明显差异时才返回
        if best_ratio > 0:
            return best_color
        return None

    def get_colors(self, h1, h2) -> Set[Tuple[str, int]]:
        """返回 (h1, h2) 上的颜色集。"""
        if not self.has_edge(h1, h2):
            return set()
        return set(self[h1][h2]['colors'])

    def edge_count(self) -> int:
        return self.number_of_edges()

    def node_count(self) -> int:
        return self.number_of_nodes()

    def all_hogs(self) -> Set:
        """所有 HOG 节点。"""
        return {n for n in self.nodes if not isinstance(n, tuple)
                or len(n) != 2 or n[1] not in ('L', 'R')}

    def log_events_summary(self, phase_label: str, event_types: set = None):
        """统一的事件 log 格式。按孩子分组输出。"""
        from collections import defaultdict
        events = self.events
        if event_types:
            events = [e for e in events if e.event_type in event_types]
        if not events:
            return
        child_events = defaultdict(list)
        for e in events:
            parts = e.branch.split('-')
            cid = parts[-1] if parts else 'unknown'
            child_events[cid].append(e)
        for cid in sorted(child_events.keys()):
            evts = child_events[cid]
            type_counts = defaultdict(int)
            type_lens = defaultdict(list)
            for e in evts:
                type_counts[e.event_type] += 1
                type_lens[e.event_type].append(len(e.genes_involved))
            parts = []
            for etype in sorted(type_counts.keys()):
                count = type_counts[etype]
                lens = type_lens[etype]
                if len(lens) == 1:
                    parts.append(f"{etype}={count} (len {lens[0]})")
                else:
                    parts.append(f"{etype}={count} (len {min(lens)}-{max(lens)})")
            logger.info("  [%s] %s events: %s", phase_label, cid, ", ".join(parts))

    def shared_edges(self) -> List[Tuple]:
        """返回多于一个颜色的边 → 祖先共享邻接。"""
        result = []
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) > 1:
                result.append((h1, h2))
        return result

    def unique_edges(self) -> List[Tuple]:
        """返回只有一个颜色的边 → 可能为衍生边。"""
        result = []
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) == 1:
                result.append((h1, h2))
        return result

    def edges_by_color(self, color: Tuple[str, int]) -> List[Tuple]:
        """返回所有带有指定颜色的边。"""
        result = []
        for h1, h2, data in self.edges(data=True):
            if color in data['colors']:
                result.append((h1, h2))
        return result

    def children(self) -> Set[str]:
        """所有参与颜色的 child_id 集合。"""
        children = set()
        for _, _, data in self.edges(data=True):
            for child_id, _ in data['colors']:
                children.add(child_id)
        return children

    # ==================== 方向调和 ====================

    def harmonize_directions(self):
        """方向调和：选一个孩子作参考，翻转方向相反的孩子。

        按共享边内容匹配（不依赖 chrom_idx）：
        - 找两个孩子共有的边
        - 统计方向一致/相反的比例
        - 多数边方向相反 → 翻转该孩子的所有边

        整条染色体反向 → 翻转方向（不是倒位事件）。
        部分段反向 → 保留方向冲突（倒位事件，后续检测）。
        """
        children = list(self.children())
        if len(children) < 2:
            return

        ref_cid = children[0]

        # 收集参考孩子的边方向
        # ref_dir: (h1, h2) -> direction
        ref_dir = {}
        for h1, h2, data in self.edges(data=True):
            for cid, ci, d in data.get('directions', set()):
                if cid == ref_cid:
                    ref_dir[(h1, h2)] = d
                    break

        # 对每个其他孩子，比较共享边方向
        flipped = []  # (cid, agree, disagree, shared, chroms, edges)
        for cid in children[1:]:
            child_dir = {}
            for h1, h2, data in self.edges(data=True):
                for c, ci, d in data.get('directions', set()):
                    if c == cid:
                        child_dir[(h1, h2)] = d
                        break

            # 找共享边并比较方向
            shared_edges = set(ref_dir.keys()) & set(child_dir.keys())
            if not shared_edges:
                continue

            agree = 0
            disagree = 0
            for edge in shared_edges:
                if ref_dir[edge] == child_dir[edge]:
                    agree += 1
                else:
                    disagree += 1

            if disagree > agree:
                n_chroms, n_edges = self._flip_all_directions(cid)
                flipped.append((cid, agree, disagree, len(shared_edges), n_chroms, n_edges))

        if flipped:
            parts = ", ".join(f"{cid}({str_agg}/{str_dis}/{str_shr}→{nc}c/{ne}e)"
                              for cid, str_agg, str_dis, str_shr, nc, ne in flipped)
            logger.info("  [harmonize] %d children flipped: %s", len(flipped), parts)

    def _flip_all_directions(self, child_id: str):
        """翻转指定孩子的所有边方向。返回 (chroms_flipped, edges_flipped)。"""
        edges_flipped = 0
        chroms = set()
        for h1, h2, data in self.edges(data=True):
            directions = data.get('directions', set())
            new_directions = set()
            for cid, ci, d in directions:
                if cid == child_id:
                    new_directions.add((cid, ci, -d))
                    chroms.add(ci)
                else:
                    new_directions.add((cid, ci, d))
            if new_directions != directions:
                data['directions'] = new_directions
                edges_flipped += 1
        return len(chroms), edges_flipped

    def remove_edge_color(self, h1, h2, color: Tuple[str, int]):
        """移除边上的一个颜色。如果该边没有其他颜色，移除整条边。"""
        if not self.has_edge(h1, h2):
            return
        colors = self[h1][h2]['colors']
        colors.discard(color)
        if not colors:
            self.remove_edge(h1, h2)

    def remove_edges_with_color(self, color: Tuple[str, int]):
        """移除所有带有指定颜色的边。"""
        to_remove = list(self.edges_by_color(color))
        for h1, h2 in to_remove:
            self.remove_edge_color(h1, h2, color)

    # ==================== Indel/Loss 检测 ====================

    def find_spanning_edges(self) -> List[Tuple]:
        """找跨越边: unique_edge (a, b) 的端点在另一孩子的路径中连通。

        对每条 unique_edges，检查 a 和 b 是否在同一孩子的另一条路径上:
        - 获取该 unique_edge 的颜色
        - 对于其他孩子的图，检查 a 和 b 之间是否有路径（且路径长度 > 1）

        Returns: [(h1, h2, child_id, spanned_hogs), ...]
        """
        spanning = []
        for h1, h2 in self.unique_edges():
            colors = self.get_colors(h1, h2)
            if not colors:
                continue
            child_id, _ = next(iter(colors))

            # 对每个其他孩子，检查是否有跨越关系
            for other_id in self.children():
                if other_id == child_id:
                    continue
                other_hogs = self._child_hogs.get(other_id, set())
                if h1 not in other_hogs or h2 not in other_hogs:
                    continue

                # 检查 h1 和 h2 是否在 other 的图中通过中间节点连通
                path = self._shortest_path_through_hogs(h1, h2, other_id)
                if path and len(path) > 2:  # 中间至少有一个 HOG
                    spanned = path[1:-1]  # 跳过两端
                    spanning.append((h1, h2, child_id, spanned))
                    break
        return spanning

    def _shortest_path_through_hogs(self, a, b, child_id: str) -> Optional[List]:
        """在指定孩子的子图中找 a→b 的最短路径（仅用该孩子的边）。"""
        G_sub = nx.DiGraph()
        for h1, h2, data in self.edges(data=True):
            if any(cid == child_id for cid, _ in data['colors']):
                G_sub.add_edge(h1, h2)
        try:
            path = nx.shortest_path(G_sub, a, b)
            return path
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            return None

    # ==================== 环检测 + 事件分类 ====================

    def find_cycles(self) -> List[List]:
        """nx.cycle_basis 包装。"""
        try:
            return nx.cycle_basis(self)
        except Exception:
            return []

    def classify_cycle(self, cycle: List) -> Tuple[Optional[str], List, Dict]:
        """分析环的颜色模式，判断事件类型。

        方向感知: 先检查环中是否有方向冲突的边对。
        - 有方向冲突 → 真实结构事件 (inversion/RT)
        - 无方向冲突 → indel shortcut 假环 (不应存在，标记为 gene_indel)

        Returns:
            (event_type, conflict_edges, info_dict)
            event_type = None 表示无法分类
            conflict_edges = 需要移除的边列表 [(h1, h2), ...]
            info_dict = 用于调试的额外信息
        """
        n = len(cycle)
        if n < 3:
            return None, [], {}

        # 构建环的边序列
        edges = [(cycle[i], cycle[(i + 1) % n]) for i in range(n)]
        edge_colors = [self.get_colors(u, v) for u, v in edges]

        info = {
            'cycle': cycle,
            'edges': edges,
            'edge_colors': edge_colors,
            'n_unique': sum(1 for c in edge_colors if len(c) == 1),
            'n_shared': sum(1 for c in edge_colors if len(c) > 1),
        }

        # == 方向冲突检查: 区分真实环 (inversion/RT) 和 indel 假环 ==
        has_dir_conflict = False
        for (u, v) in edges:
            if self.edge_has_direction_conflict(u, v):
                has_dir_conflict = True
                break
        info['has_direction_conflict'] = has_dir_conflict

        if not has_dir_conflict:
            # 无方向冲突 → indel shortcut 假环
            # 这些边不应该形成环，标记为 gene_indel 让 resolve_indels 处理
            info['resolution'] = 'no_direction_conflict_indel'
            return 'gene_indel', list(edges), info

        # 统计每种颜色的出现次数
        color_counts = Counter()
        for colors in edge_colors:
            for c in colors:
                color_counts[c] += 1

        info['color_counts'] = dict(color_counts)

        # 如果颜色数 = 0 (无颜色)，无法分类
        if not color_counts:
            return None, [], info

        # == 情况 A: 4-HOG 环，颜色交替 (A1, B1, A1, B1) ==
        if n >= 4:
            # 检查是否交替: edges[i] 和 edges[i+2] 应有相同的颜色集
            alternating = True
            for i in range(n):
                j = (i + 2) % n
                if edge_colors[i] != edge_colors[j]:
                    alternating = False
                    break

            if alternating and n >= 4:
                # 交替模式 → inversion 或 RT/URT
                # 检查是否涉及不同染色体
                chroms = set()
                for colors in edge_colors:
                    for cid, chrom_idx in colors:
                        chroms.add((cid, chrom_idx))
                # 获取边缘的颜色集
                col_set = set()
                for colors in edge_colors:
                    col_set.update(colors)

                if len(col_set) >= 3:
                    # 3+ 种颜色 → 跨染色体 → RT/URT
                    # 冲突边: 最少出现的颜色的边
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    conflict_edges = []
                    for (h1, h2), colors in zip(edges, edge_colors):
                        if rare_colors & colors:
                            conflict_edges.append((h1, h2))
                    info['conflict_colors'] = rare_colors
                    # 检查是否含端粒 → URT
                    has_tel = any(
                        isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R')
                        for n in cycle
                    )
                    etype = 'unbalanced_reciprocal_translocation' if has_tel else 'reciprocal_translocation'
                    return etype, conflict_edges, info
                else:
                    # 2 种颜色 → inversion
                    # 端粒驱动：保留端粒邻接的颜色，移除另一种
                    tel_color = self.telomere_preserving_color(edges, edge_colors)

                    if tel_color is not None:
                        # 端粒信号明确：移除不保留端粒的颜色
                        remove_colors = {c for c in color_counts if c != tel_color}
                        info['resolution'] = 'telomere-driven'
                        info['kept_color'] = tel_color
                    else:
                        # 无端粒信号：用外类群或回退到出现次数
                        min_count = min(color_counts.values())
                        remove_colors = {c for c, cnt in color_counts.items()
                                         if cnt == min_count}
                        if len(remove_colors) > 1:
                            remove_colors = {next(iter(remove_colors))}
                        info['resolution'] = 'count-based fallback'

                    conflict_edges = []
                    for (h1, h2), colors in zip(edges, edge_colors):
                        if remove_colors & colors:
                            conflict_edges.append((h1, h2))
                    info['conflict_colors'] = remove_colors

                    has_tel = any(
                        isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R')
                        for n in cycle
                    )
                    etype = 'telomere_inversion' if has_tel else 'inversion'
                    return etype, conflict_edges, info

        # == 简单情况: 3-HOG 环 (不可能在染色体图中, 属于错误) ==
        if n == 3:
            return 'gene_indel', list(edges), info

        # 无法分类的环 → 优先用端粒信号，否则移除最少出现的颜色
        tel_color = self.telomere_preserving_color(edges, edge_colors)
        if tel_color is not None:
            remove_colors = {c for c in color_counts if c != tel_color}
            info['resolution'] = 'telomere-driven fallback'
        else:
            min_count = min(color_counts.values())
            remove_colors = {c for c, cnt in color_counts.items() if cnt == min_count}
            info['resolution'] = 'count-based fallback'

        conflict_edges = []
        for (h1, h2), colors in zip(edges, edge_colors):
            if remove_colors & colors:
                conflict_edges.append((h1, h2))
        info['fallback'] = True
        info['fallback_reason'] = 'unclassified_cycle'
        return 'unclassified', conflict_edges, info

    # ==================== 结构重排 resolved ====================

    def resolve_structural_events(self, outgroup_adjacency: Set = None):
        """循环: find_cycles → classify_cycle → 移除外类群确认的衍生方边。

        Outgroup voting 确定极性:
        对环中每条边 (hi, hj), 映射到 parent HOG level 查外类群邻接集。
        - 外类群有该邻接 → 祖先状态 → 保存该边的孩子是祖先方
        - 外类群无该邻接 → 衍生状态 → 有该边的孩子是衍生方
        - 衍生方 = 事件发生方 → 移除该孩子的边

        无外类群 (root): 沿用最少颜色出现次数的启发式。
        """
        iteration = 0
        max_iterations = 100  # 防止无限循环
        while iteration < max_iterations:
            iteration += 1
            cycles = self.find_cycles()
            if not cycles:
                break
            logger.debug("  [structural] iteration %d: %d cycles found",
                         iteration, len(cycles))
            resolved_this_round = 0
            for cycle in cycles:
                etype, conflict_edges, info = self.classify_cycle(cycle)
                if etype is None or not conflict_edges:
                    continue

                # 构建环的边序列（用于外类群投票和回退模式）
                n_cycle = len(cycle)
                cycle_edges = [(cycle[i], cycle[(i + 1) % n_cycle])
                               for i in range(n_cycle)]

                # 用外类群投票确定衍生方（事件发生分支）
                if outgroup_adjacency is not None and etype in (
                    'inversion', 'telomere_inversion',
                    'reciprocal_translocation',
                    'unbalanced_reciprocal_translocation',
                ):
                    cycle_colors = [self.get_colors(u, v) for u, v in cycle_edges]

                    # 对每条边，映射到 parent HOG level 查外类群
                    child_derived_count = defaultdict(int)
                    child_ancestral_count = defaultdict(int)
                    for (u, v), colors in zip(cycle_edges, cycle_colors):
                        u_parent = getattr(u, 'parent',
                                           u.hog_id if hasattr(u, 'hog_id') else str(u))
                        v_parent = getattr(v, 'parent',
                                           v.hog_id if hasattr(v, 'hog_id') else str(v))
                        key = (u_parent, v_parent) if u_parent < v_parent else (v_parent, u_parent)
                        is_ancestral = key in outgroup_adjacency
                        for cid, _ in colors:
                            if is_ancestral:
                                child_ancestral_count[cid] += 1
                            else:
                                child_derived_count[cid] += 1

                    # 出现衍生边的孩子 = 事件发生方
                    if child_derived_count:
                        derived_child_id = max(child_derived_count,
                                               key=child_derived_count.get)
                    else:
                        derived_child_id = None

                    logger.debug("  [structural] cycle %s: og-derived=%s, "
                                 "ancestral=%s, derived=%s",
                                 cycle, derived_child_id,
                                 dict(child_ancestral_count),
                                 dict(child_derived_count))

                    # 只移除衍生方的边
                    if derived_child_id:
                        for h1, h2 in conflict_edges:
                            colors = self.get_colors(h1, h2)
                            for cid, chrom in list(colors):
                                if cid == derived_child_id:
                                    self.remove_edge_color(h1, h2, (cid, chrom))
                    else:
                        # 所有边都是祖先状态？回退到最少颜色出现次数的启发式
                        for h1, h2 in conflict_edges:
                            self._remove_rare_color_edge(h1, h2, cycle, cycle_edges)
                else:
                    # 无外类群 (root) 或 unclassified 环 → 启发式
                    for h1, h2 in conflict_edges:
                        self._remove_rare_color_edge(h1, h2, cycle, cycle_edges)

                # 记录事件
                _dcid = derived_child_id if outgroup_adjacency and child_derived_count else None
                event = TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{_dcid}" if _dcid else self.hog_level,
                    genes_involved=list(cycle),
                    desc=f"{etype}: {len(cycle)} HOG cycle resolved (iter {iteration})",
                    support=1,
                )
                self.events.append(event)
                resolved_this_round += 1

            if resolved_this_round == 0:
                logger.warning("  [structural] %d cycles found but none resolved",
                               len(cycles))
                break

        if iteration >= max_iterations:
            logger.warning("  [structural] max iterations reached (%d), force-breaking cycles",
                           max_iterations)
            # 强制打破剩余环：对每个环移除最少 child 的所有边
            cycles = self.find_cycles()
            for cycle in cycles:
                etype, conflict_edges, info = self.classify_cycle(cycle)
                for h1, h2 in conflict_edges:
                    colors = self.get_colors(h1, h2)
                    if len(colors) == 1:
                        self.remove_edge_color(h1, h2, next(iter(colors)))
                    elif colors:
                        self.remove_edge_color(h1, h2, next(iter(colors)))
                self.events.append(TAKREvent(
                    event_type='unclassified',
                    branch=self.hog_level,
                    genes_involved=list(cycle),
                    desc=f"force-break: {len(cycle)} HOGs",
                    support=1,
                ))

    # ==================== 边移除辅助 ====================

    def _remove_rare_color_edge(self, h1, h2, cycle, cycle_edges):
        """移除冲突边上出现次数最少的颜色（回退启发式）。"""
        colors = self.get_colors(h1, h2)
        if len(colors) == 1:
            self.remove_edge_color(h1, h2, next(iter(colors)))
        else:
            # 统计各 (cid, chrom) 在环中的出现次数
            child_count = {}
            for (u, v) in cycle_edges:
                for cid, chrom in self.get_colors(u, v):
                    child_count[(cid, chrom)] = child_count.get((cid, chrom), 0) + 1
            if child_count:
                min_child = min(child_count, key=child_count.get)
                if min_child in colors:
                    self.remove_edge_color(h1, h2, min_child)
                else:
                    self.remove_edge_color(h1, h2, next(iter(colors)))
            elif colors:
                self.remove_edge_color(h1, h2, next(iter(colors)))

    # ==================== 路径覆盖 ====================

    def path_cover(self) -> List[List]:
        """端粒约束路径覆盖。

        生物学约束：每条祖先染色体有且仅有 2 个端粒（两端各一个）。
        端粒位置通过多个孩子的共识确定。

        算法：
        1. 找共识端粒 HOG（≥2 个孩子中都与端粒相邻）
        2. 从共识端粒出发，沿邻接行走
        3. 行走时优先选择端粒邻接度高的邻居
        4. 到达另一个端粒时停止 → 一条完整染色体
        5. 处理剩余节点（无端粒的连通分量）

        Returns:
            [[hog1, hog2, ...], ...] — 每条染色体一条路径
        """
        paths = []
        visited = set()
        hog_set = self.all_hogs()  # cache — avoid repeated O(n) calls

        # 共识端粒：在任一孩子中与端粒相邻的 HOG
        # (min_children=1 更宽松，覆盖 rearrangement 后的端粒)
        cons_tels = self.consensus_telomeres(min_children=1)
        # 备用：所有孩子的端粒并集
        all_tels = self.child_telomere_set()
        # 度=1 的节点（图结构上的端点）
        degree_ends = {n for n in self.nodes()
                       if self.in_degree(n) == 1 and n in hog_set}

        # 选择端粒锚点集：共识端粒 > 所有端粒 > 度=1
        # 注意：度=1 包含 bridge 移除后的假端点，优先级最低
        if cons_tels:
            anchor_tels = cons_tels
        elif all_tels:
            anchor_tels = all_tels
        else:
            anchor_tels = degree_ends

        # 从每个锚点出发行走
        for start in sorted(anchor_tels, key=str):
            if start in visited or not self.has_node(start):
                continue

            path = self._walk_telomere_path(start, visited, anchor_tels)
            if path and len(path) >= 1:
                paths.append(path)
                visited.update(path)

        # 处理剩余未访问节点（无端粒的连通分量 + 孤立节点）
        all_nodes = self.all_hogs()
        unvisited = all_nodes - visited
        if unvisited:
            G_rem = self.subgraph(unvisited)
            for comp_nodes in nx.weakly_connected_components(G_rem):
                comp = list(comp_nodes)
                if len(comp) >= 2:
                    start = comp[0]
                    path = self._walk_simple(start, visited)
                    if path:
                        paths.append(path)
                        visited.update(path)
                elif len(comp) == 1:
                    # 孤立节点：不单独成染色体（noise），跳过
                    visited.update(comp)

        # 孤儿回收：未覆盖的 HOGs 插入到最近的已覆盖路径
        # 如果无法插入，作为独立路径
        all_nodes = self.all_hogs()
        unvisited = all_nodes - visited
        if unvisited:
            logger.debug("  [path_cover] %d uncovered HOGs, attempting recovery",
                         len(unvisited))
            recovered = 0
            for hog in list(unvisited):
                if not self.has_node(hog):
                    continue
                neighbors = [n for n in self._neighbors(hog)
                             if n in visited and n in all_nodes]
                if not neighbors:
                    continue
                target_path = None
                target_idx = None
                for nb in neighbors:
                    for pi, p in enumerate(paths):
                        if nb in p:
                            idx = p.index(nb)
                            target_path = pi
                            target_idx = idx
                            break
                    if target_path is not None:
                        break
                if target_path is not None:
                    paths[target_path].insert(target_idx + 1, hog)
                    visited.add(hog)
                    recovered += 1
            if recovered:
                logger.debug("  [path_cover] recovered %d orphan HOGs", recovered)

            # 仍有未覆盖的 HOGs → 断开组件作为独立路径
            still_unvisited = all_nodes - visited
            if still_unvisited:
                G_rem = self.subgraph(still_unvisited)
                for comp_nodes in nx.weakly_connected_components(G_rem):
                    comp = list(comp_nodes)
                    if len(comp) >= 2:
                        start = comp[0]
                        path = self._walk_simple(start, visited)
                        if path:
                            paths.append(path)
                            visited.update(path)
                    elif len(comp) == 1:
                        # 孤立节点也作为路径（保留信息）
                        paths.append(list(comp))
                        visited.update(comp)

        return paths

    def _walk_telomere_path(self, start, visited: Set,
                            anchor_tels: Set) -> Optional[List]:
        """从端粒 HOG 出发，走到另一个端粒 HOG。

        行走策略：
        - 度=1 的邻居优先（端点）
        - 端粒锚点 HOG 次之
        - 遇到另一个锚点时停止
        """
        path = [start]
        curr = start
        visited_local = visited | {start}
        hog_set = self.all_hogs()  # cache

        while True:
            neighbors = [n for n in self._neighbors(curr)
                         if n not in visited_local and n in hog_set]
            if not neighbors:
                break

            if len(neighbors) == 1:
                nxt = neighbors[0]
            else:
                # 多个邻居时，按优先级排序
                def neighbor_priority(n):
                    if self._degree(n) == 1:
                        return 0
                    if n in anchor_tels:
                        return 1
                    return 2
                neighbors.sort(key=neighbor_priority)
                nxt = neighbors[0]

            visited_local.add(nxt)
            path.append(nxt)
            curr = nxt

            # 到达另一个锚点 → 停止
            if curr != start and curr in anchor_tels:
                break

        return path if len(path) >= 1 else None

    def _walk_simple(self, start, visited: Set) -> Optional[List]:
        """简单行走（无端粒约束），用于处理剩余连通分量。"""
        path = [start]
        curr = start
        visited_local = visited | {start}
        hog_set = self.all_hogs()  # cache

        while True:
            neighbors = [n for n in self._neighbors(curr)
                         if n not in visited_local and n in hog_set]
            if not neighbors:
                break
            nxt = neighbors[0]
            visited_local.add(nxt)
            path.append(nxt)
            curr = nxt

        return path if len(path) >= 1 else None

    # ==================== 转换 ====================

    def to_ancestral_graph(self):
        """转换为 AncestralAdjacencyGraph。

        复用 resolve_all_events 中缓存的 paths，避免重复计算。
        """
        from .AK import AncestralAdjacencyGraph

        result = AncestralAdjacencyGraph(node_id=self.hog_level)

        # 复用缓存的 paths
        paths = getattr(self, '_cached_paths', None) or self.path_cover()
        if not paths:
            logger.warning("  [colored] no paths found for %s", self.hog_level)
            return result

        for chrom_idx, path in enumerate(paths):
            # 添加 HOG 节点
            for n in path:
                if n not in result.graph:
                    result.graph.add_node(n)
                    result.gene_nodes.add(n)

            # 添加有向边（基于路径顺序）
            for i in range(len(path) - 1):
                result.graph.add_edge(path[i], path[i + 1])

            # 添加端粒
            chrom_name = f"{self.hog_level}_chrom_{chrom_idx}"
            left_tel = (chrom_name, 'L')
            right_tel = (chrom_name, 'R')
            result.graph.add_node(left_tel, telomere=True)
            result.graph.add_node(right_tel, telomere=True)
            result.telomeres.add(left_tel)
            result.telomeres.add(right_tel)
            result.graph.add_edge(left_tel, path[0])
            result.graph.add_edge(path[-1], right_tel)

        result.events = self.events
        return result

    def check_block_degrees(self) -> List:
        """检查 block 级有向图：非端粒 block 必须有 in>0 且 out>0。

        Returns: [(bid, in_deg, out_deg, chr_name, size), ...] 不合格的 block 列表。
        """
        bg = getattr(self, '_block_graph', None)
        bg = getattr(self, '_block_graph', None)
        if bg is None or bg.number_of_nodes() == 0:
            return []
        broken = []
        for bid in bg.nodes():
            hogs = self._blocks.get(bid, [])
            # 端粒判定：单 HOG block 且 HOG 带 telomere 属性
            is_tel = (len(hogs) == 1 and
                      self.has_node(hogs[0]) and
                      bool(self.nodes[hogs[0]].get('telomere')))
            if is_tel:
                continue
            in_d = bg.in_degree(bid)
            out_d = bg.out_degree(bid)
            if in_d == 0 or out_d == 0:
                ch = '?'
                for h in hogs[:1]:
                    srcs = self.nodes[h].get('sources', set())
                    if srcs:
                        ch = next(iter(srcs))[1]
                        break
                broken.append((bid, in_d, out_d, ch, len(hogs)))
        return broken

    def _ensure_blocks(self):
        """确保 blocks 和 block_graph 已构建。幂等。"""
        if not hasattr(self, '_blocks') or not self._blocks:
            self._build_synteny_blocks()
        if not hasattr(self, '_block_graph') or not self._block_graph:
            self._compress_to_block_level()

    def to_gfa(self, fout, min_pair_nodes: int = 50, use_blocks: bool = True,
               color_by: str = 'first'):
        """输出 GFA 格式。use_blocks=False 输出 HOG 级，否则 block 级。

        着色规则:
          - 端粒节点 → 红色 (#FF0000)，优先级最高
          - 非端粒节点 → color_by='first': 按第一个物种染色体着色（默认）
                          color_by='all': 按所有物种|染色体组合着色
          - 每个节点输出 LB:Z:label 标明来源

        GFA 格式:
          S  node  *  CL:Z:#RRGGBB  LB:Z:label  LN:i:N
          L  n1  +  n2  +  0M  CO:Z:c1:0,c2:1,...
        """
        from collections import defaultdict

        # ── 确保 blocks 已构建 ──
        if use_blocks:
            self._ensure_blocks()
        if use_blocks:
            block_graph = getattr(self, '_block_graph', None)
            graph = block_graph if (block_graph is not None
                                    and block_graph.number_of_nodes() > 0) else self
        else:
            graph = self
            block_graph = None
        _blocks = getattr(self, '_blocks', {})

        # ── 收集节点 label 和计数 ──
        label_counts = defaultdict(int)
        node_label = {}
        node_sources = {}   # LB 标签：完整来源信息

        for node in graph.nodes:
            sources = graph.nodes[node].get('sources', set())
            is_tel = bool(graph.nodes[node].get('telomere'))

            full_src = '+'.join(sorted(f'{cid}|{ch}' for cid, ch in sources)) if sources else ''
            node_sources[node] = full_src

            if sources:
                if color_by == 'first':
                    first_sp = sorted(sources)[0][0]
                    chrom_label = '+'.join(sorted(ch for cid, ch in sources if cid == first_sp))
                    label = chrom_label if chrom_label else full_src
                else:
                    label = full_src
            else:
                label = ''
            node_label[node] = label
            if label and not is_tel:
                n_len = len(_blocks.get(node, [])) if use_blocks else 1
                label_counts[label] += n_len

        # ── 色板 ──
        significant = {l for l, c in label_counts.items() if c >= min_pair_nodes}
        label_colors = _assign_pair_colors(sorted(significant)) if significant else {}

        # ── Segment 行 ──
        for node in sorted(graph.nodes, key=str):
            label = node_label[node]
            if use_blocks:
                # 端粒 block: 必须是单 HOG 且带 telomere 属性
                is_tel = (len(_blocks.get(node, [])) == 1
                          and bool(graph.nodes[node].get('telomere')))
            else:
                is_tel = bool(self.nodes[node].get('telomere'))

            if is_tel:
                color = '#FF0000'
            elif label and label in label_colors:
                color = label_colors[label]
            else:
                color = '#808080'

            line = ['S', node, '*', f'CL:Z:{color}']
            src = node_sources.get(node, '')
            lb = f'{src}+telo' if (is_tel and src) else ('telo' if is_tel else src)
            if lb:
                line.append(f'LB:Z:{lb}')
            # Length: 端粒固定 50（显眼），否则 block 大小
            n_len = 50 if is_tel else (len(_blocks.get(node, [])) if use_blocks else 1)
            line.append(f'LN:i:{n_len}')
            fout.write('\t'.join(map(str, line)) + '\n')

        # ── A 行：block → HOG 映射（标准 GFA alignment 格式）──
        if use_blocks:
            for bid in sorted(_blocks.keys(), key=lambda x: int(x.split('_')[1]) if '_' in x else 0):
                for pos, hog in enumerate(_blocks[bid]):
                    fout.write(f'A\t{bid}\t{pos}\t+\t{hog}\t0\t1\n')

        # ── Link 行 ──
        for n1, n2, data in sorted(graph.edges(data=True),
                                    key=lambda x: (str(x[0]), str(x[1]))):
            colors = data.get('colors', set())
            color_tags = ','.join(f'{c[0]}:{c[1]}' for c in sorted(colors))
            line = ['L', n1, '+', n2, '+', '0M']
            if color_tags:
                line.append(f'CO:Z:{color_tags}')
            fout.write('\t'.join(map(str, line)) + '\n')

    # ==================== 一键执行 ====================

    def collapse_wgd(self, ploidy: int):
        """ColoredGraph-based WGD collapse: post-WGD → pre-WGD.

        post-WGD 图有 ploidy×pre 条染色体，它们是天然分开的（无跨染色体边）。
        WGD collapse = 染色体配对：用 HOG 内容相似度找亚基因组同位染色体。

        Pipeline:
        1. Path_cover → 获取 N 条染色体
        2. 染色体配对: 每对染色体的 HOG Jaccard 相似度 → 找最佳匹配
        3. 每对合并为一条 pre-WGD 染色体（HOG 去重）

        Args:
            ploidy: 倍性 (2=四倍体, 3=六倍体, ...)

        Returns:
            ColoredGraph — pre-WGD 图
        """
        # Step 1: 获取 post-WGD 染色体
        paths = self.path_cover()
        if not paths:
            logger.warning("  [collapse_wgd] no paths found, returning original graph")
            return self

        n_chrom = len(paths)
        logger.debug("  [collapse_wgd] %s: %d post-WGD chromosomes (ploidy=%d)",
                      self.hog_level, n_chrom, ploidy)

        if n_chrom < ploidy:
            logger.warning("  [collapse_wgd] fewer chroms (%d) than ploidy (%d), skipping",
                           n_chrom, ploidy)
            return self

        # Step 2: 构建染色体 × HOG 矩阵
        chrom_hogs = []  # [{hog_id}, ...] per chromosome
        for path in paths:
            chrom_hogs.append(set(path))

        # 计算所有染色体对的 Jaccard 相似度
        n = len(chrom_hogs)
        import itertools
        pairs = []
        for i, j in itertools.combinations(range(n), 2):
            inter = len(chrom_hogs[i] & chrom_hogs[j])
            union = len(chrom_hogs[i] | chrom_hogs[j])
            jaccard = inter / union if union > 0 else 0.0
            if inter > 0:
                pairs.append((jaccard, inter, i, j))

        # Step 3: 贪心匹配 — 从相似度最高的对开始
        pairs.sort(reverse=True)
        paired = set()
        pre_chroms = []  # [(path_i, path_j), ...]

        for jaccard, inter, i, j in pairs:
            if i in paired or j in paired:
                continue
            pre_chroms.append((paths[i], paths[j]))
            paired.add(i)
            paired.add(j)

        # 未配对的染色体作为单条保留
        for i in range(n):
            if i not in paired:
                pre_chroms.append((paths[i],))

        n_pre = len(pre_chroms)
        logger.debug("  [collapse_wgd] %s: %d post-WGD -> %d pre-WGD chroms (paired=%d)",
                      self.hog_level, n_chrom, n_pre, len(pre_chroms))

        # Step 4: 构建 pre-WGD 图（每对染色体合并）
        G2 = ColoredGraph(hog_level=f"{self.hog_level}_preWGD")

        # 记录事件
        paired_count = 0
        for chrom_pair in pre_chroms:
            if len(chrom_pair) == 2:
                # 一对亚基因组同位染色体 → 合并（HOG 去重）
                path_a, path_b = chrom_pair
                merged = _merge_chromosome_paths(path_a, path_b)
                for i in range(len(merged) - 1):
                    G2.add_edge(merged[i], merged[i + 1],
                                f"pre_chr{paired_count}", paired_count)
                paired_count += 1
                # 记录配对事件
                shared_hogs = set(path_a) & set(path_b)
                if shared_hogs:
                    G2.events.append(TAKREvent(
                        event_type='wgd_collapse',
                        branch=f"{self.hog_level}_preWGD-{self.hog_level}",
                        genes_involved=list(shared_hogs),
                        desc=f"wgd_collapse: paired chroms ({len(path_a)}+{len(path_b)} HOGs, {len(shared_hogs)} shared)",
                        support=ploidy,
                    ))
            else:
                # 未配对的单条染色体
                path = chrom_pair[0]
                for i in range(len(path) - 1):
                    G2.add_edge(path[i], path[i + 1],
                                f"pre_chr{paired_count}", paired_count)
                paired_count += 1

        return G2

    # ==================== 共线性块压缩 ====================

    # ==================== 共线性块压缩 ====================

    @staticmethod
    def _extract_unitigs(G, min_size=2, skip_nodes=None):
        """有向图 unitig：每个孩子 indegree=1 且 outdegree=1 的连续节点。

        skip_nodes 中的节点为 stop——不进入 unitig。
        相邻节点物种来源不一致也断——单孩子/多孩子统一处理。
        双向边按单孩子视角判断 linearity，两个孩子方向相反也各自算 linear。
        """
        skip = skip_nodes or set()

        def _species_set(node):
            return set(cid for cid, _ in G.nodes[node].get('sources', set()))

        def _linear(node, direction):
            """所有孩子是否同意同一个 direction 方向上的邻居。

            direction='succ': fwd walk 可继续
            direction='pred': bwd walk 可继续

            定向不一致时回退到邻居对 frozenset 检查——
            允许双向边（孩子 A 的 succ = 孩子 B 的 pred，pair 相同）。
            """
            if G.nodes[node].get('telomere') or node in skip:
                return False
            deg = ColoredGraph._per_child_degree(G, node)
            if not deg:
                return False
            for (inn, out) in deg.values():
                if inn != 1 or out != 1:
                    return False
            sources = G.nodes[node].get('sources', set())
            children = set(c for c, _ in sources)
            if len(children) <= 1:
                return True
            # 定向一致性检查
            neighbors = list(G.successors(node) if direction == 'succ'
                             else G.predecessors(node))
            agreed = None
            for cid in children:
                child_nbrs = [n for n in neighbors
                              if any(c == cid for c, _ in
                                     (G[node][n] if direction == 'succ'
                                      else G[n][node]).get('colors', set()))]
                if len(child_nbrs) != 1:
                    return False
                if agreed is None:
                    agreed = child_nbrs[0]
                elif agreed != child_nbrs[0]:
                    # 定向不一致 → 回退到 pair 检查（双向边）
                    break
            else:
                return True  # 所有孩子定向一致
            # 回退：邻居对 frozenset 一致（双向边兼容）
            succs = list(G.successors(node))
            preds = list(G.predecessors(node))
            pair0 = None
            for cid in children:
                s = next((n for n in succs
                          if any(c == cid for c, _ in G[node][n].get('colors', set()))), None)
                p = next((n for n in preds
                          if any(c == cid for c, _ in G[n][node].get('colors', set()))), None)
                if s is None or p is None:
                    return False
                pair = frozenset([s, p])
                if pair0 is None:
                    pair0 = pair
                elif pair != pair0:
                    return False
            return True

        unitigs = []
        visited = set()
        for node in G.nodes():
            if node in visited:
                continue

            # 从当前节点出发，向前走后继，向回走前驱
            fwd = []
            curr = node
            while True:
                if G.nodes[curr].get('telomere') or curr in skip:
                    break
                deg = ColoredGraph._per_child_degree(G, curr)
                if not deg or any(inn != 1 or out != 1
                                  for (inn, out) in deg.values()):
                    break
                fwd.append(curr)
                visited.add(curr)
                if not _linear(curr, 'succ'):
                    break
                nxt = ColoredGraph._per_child_neighbor(G, curr, 'succ')
                if nxt is None:
                    break
                # 后继已访问或物种不一致 → 断（curr 已在 unitig 内）
                if nxt in visited or _species_set(curr) != _species_set(nxt):
                    break
                curr = nxt

            bwd = []
            curr = node
            while True:
                nxt = ColoredGraph._per_child_neighbor(G, curr, 'pred')
                if nxt is None:
                    break
                # 前驱已访问或物种不一致 → 断
                if nxt in visited or _species_set(curr) != _species_set(nxt):
                    break
                curr = nxt
                deg = ColoredGraph._per_child_degree(G, curr)
                if not deg or any(inn != 1 or out != 1
                                  for (inn, out) in deg.values()):
                    break
                if G.nodes[curr].get('telomere') or curr in skip:
                    break
                bwd.append(curr)
                visited.add(curr)
                if not _linear(curr, 'pred'):
                    break

            path = list(reversed(bwd)) + fwd
            if len(path) >= min_size:
                unitigs.append(path)
        return unitigs

    @staticmethod
    def _per_child_degree(G, node):
        """返回 {child_id: (indegree, outdegree)} 单孩子视角度数。"""
        deg = {}
        sources = G.nodes[node].get('sources', set())
        for cid in set(c for c, _ in sources):
            out = sum(1 for s in G.successors(node)
                      if any(c == cid for c, _ in G[node][s].get('colors', set())))
            inn = sum(1 for p in G.predecessors(node)
                      if any(c == cid for c, _ in G[p][node].get('colors', set())))
            deg[cid] = (inn, out)
        return deg

    @staticmethod
    def _per_child_neighbor(G, node, direction='succ'):
        """所有孩子指向同一个后继/前驱则返回该节点，否则 None。"""
        sources = G.nodes[node].get('sources', set())
        children = set(c for c, _ in sources)
        if not children:
            return None
        agreed = None
        neighbors = list(G.successors(node) if direction == 'succ'
                         else G.predecessors(node))
        for cid in children:
            child_nbrs = [n for n in neighbors
                          if any(c == cid for c, _ in
                                 (G[node][n] if direction == 'succ'
                                  else G[n][node]).get('colors', set()))]
            if len(child_nbrs) != 1:
                return None
            if agreed is None:
                agreed = child_nbrs[0]
            # 方向不一致：取第一个孩子的邻居继续（双向边进 block）
        return agreed

    def _build_synteny_blocks(self, min_block_size: int = 2):
        """构建共线性块。

        规则：
        - 端粒 HOG 永远不进入任何块（始终是 1-HOG 块）
        - 非端粒 HOG：共享边连接的连续路径 = 一个块
        - 断开点：唯一边、端粒 HOG
        """
        # 端粒 HOG 永远不压缩进 block，不论度数
        cons_tels = self.child_telomere_set()
        hog_set = self.all_hogs()
        blocks = {}
        hog_to_block = {}

        def _add_block(path):
            bid = "blk_{}".format(len(blocks))
            # 只加入尚未分配的 HOG
            deduped = [h for h in path if h not in hog_to_block]
            if not deduped:
                return None
            blocks[bid] = deduped
            for h in deduped:
                hog_to_block[h] = bid
            return bid

        # 统一调用：物种来源一致判定自动区分单/多孩子/共享段/独有段
        for path in self._extract_unitigs(self, min_block_size,
                                          skip_nodes=cons_tels):
            _add_block(path)

        # 剩余孤立 HOG → 1-HOG 块（含端粒）
        for h in hog_set:
            if h not in hog_to_block:
                _add_block([h])

        self._blocks = blocks
        self._hog_to_block = hog_to_block

        n_blk = len(blocks)
        n_multi = sum(1 for p in blocks.values() if len(p) >= 2)
        n_single = sum(1 for p in blocks.values() if len(p) == 1)
        logger.debug("  [blocks] %d blocks (%d multi-HOG, %d singleton)",
                     n_blk, n_multi, n_single)
        return blocks

    def _compress_to_block_level(self):
        """将 HOG 级图压缩为块级图（ColoredGraph）。

        直接 HOG 边 → block 边映射。
        """
        if not hasattr(self, '_blocks') or not self._blocks:
            self._build_synteny_blocks()

        block_cg = type(self)(hog_level=self.hog_level)
        for bid in self._blocks:
            block_cg.add_node(bid)
            # 全量继承第一个 HOG 的节点属性
            block_cg.nodes[bid].update(self.nodes[self._blocks[bid][0]])

        for h1, h2, data in self.edges(data=True):
            b1 = self._hog_to_block.get(h1)
            b2 = self._hog_to_block.get(h2)
            if b1 and b2:
                # 同 block 边：仅保留 HOG 级自环（h1==h2），
                # h1≠h2 的内部边已被压缩在 block 内，不映射
                if b1 == b2 and h1 != h2:
                    continue
                # 只有末端 HOG 的边才映射为块间边
                blk1 = self._blocks.get(b1, [])
                blk2 = self._blocks.get(b2, [])
                if not (blk1 and blk2
                        and (h1 == blk1[0] or h1 == blk1[-1])
                        and (h2 == blk2[0] or h2 == blk2[-1])):
                    continue
                if not block_cg.has_edge(b1, b2):
                    # 全量继承 HOG 边的所有属性
                    block_cg.add_edge(b1, b2, **data)

        self._block_graph = block_cg
        self._validate_block_compression(block_cg)

        logger.debug("  [blocks] block graph: %d nodes, %d edges",
                     block_cg.number_of_nodes(), block_cg.number_of_edges())
        return block_cg

    def _validate_block_compression(self, block_cg):
        """简单比较：hog_n = blk_n + Σ(L-1)。"""
        hog_n = self.number_of_edges()
        blk_n = block_cg.number_of_edges()
        internal = sum(len(hogs) - 1 for hogs in self._blocks.values())
        expected = blk_n + internal
        if hog_n == expected:
            logger.info("  [blocks] verified hog-blk: %d = %d + %d ✓", hog_n, blk_n, internal)
        else:
            logger.info("  [blocks] not verified : hog %d !=  blk %d + %d (diff=%d)",
                        hog_n, blk_n, internal, hog_n - expected)

    def _detect_inversions(self) -> int:
        """直接检测倒位：找方向冲突的边对。

        有符号邻接模型：同一对 HOG 在不同孩子中方向相反 → 倒位。
        不需要环检测，直接检查每条共享边的方向一致性。

        Returns: 检测到的倒位事件数
        """
        events_before = len(self.events)
        hog_set = self.all_hogs()

        for h1, h2 in list(self.edges()):
            if h1 not in hog_set or h2 not in hog_set:
                continue
            if not self.edge_has_direction_conflict(h1, h2):
                continue

            # 方向冲突 → 倒位信号
            # 确定哪个孩子是倒位方（符号与多数不同的那个）
            directions = self.get_directions(h1, h2)
            child_signs = defaultdict(set)  # child_id → set of signs
            for cid, chrom, d in directions:
                child_signs[cid].add(d)

            # 找少数方（符号与其他孩子不同的）
            all_signs = set()
            for signs in child_signs.values():
                all_signs.update(signs)

            # 每个孩子是否与其他孩子方向一致
            for cid, signs in child_signs.items():
                other_signs = set()
                for ocid, osigns in child_signs.items():
                    if ocid != cid:
                        other_signs.update(osigns)
                if signs and other_signs and not signs & other_signs:
                    # cid 的方向与其他所有孩子都不同 → cid 是倒位方
                    self.events.append(TAKREvent(
                        event_type='inversion',
                        branch=f"{self.hog_level}-{cid}",
                        genes_involved=[h1, h2],
                        desc=f"inversion: {h1}↔{h2} direction conflict in {cid}",
                        support=1,
                    ))
                    logger.debug("  [inversion] %s↔%s: child %s has opposite direction",
                                 h1, h2, cid)

        n_found = len(self.events) - events_before
        if n_found:
            logger.info("  [inversion] detected %d inversions from direction conflicts", n_found)
        return n_found

    def _find_block_cycles(self) -> List[List]:
        """在块级图上找环。"""
        if not hasattr(self, '_block_graph'):
            return []
        try:
            return nx.cycle_basis(self._block_graph.to_undirected())
        except Exception:
            return []

    def _classify_block_cycle(self, cycle: List) -> Tuple[Optional[str], List, Dict]:
        """分析块级环的颜色模式，判断事件类型。

        块级环 vs HOG 级环：
        - 环节点 = 块（block_id），而不是 HOG
        - 无需考虑块内 HOG 的顺序差异
        - 块间边的颜色 = 两个孩子对该块间邻接的贡献

        Returns:
            (event_type, conflict_block_edges, info)
        """
        n = len(cycle)
        if n < 3:
            return None, [], {}

        bg = self._block_graph

        # 构建环的边序列
        edges = [(cycle[i], cycle[(i + 1) % n]) for i in range(n)]
        edge_colors = []
        for b1, b2 in edges:
            if bg.has_edge(b1, b2):
                edge_colors.append(set(bg[b1][b2].get('colors', set())))
            else:
                edge_colors.append(set())

        info = {
            'cycle': cycle,
            'edges': edges,
            'edge_colors': edge_colors,
            'n_unique': sum(1 for c in edge_colors if len(c) == 1),
            'n_shared': sum(1 for c in edge_colors if len(c) > 1),
        }

        # 统计每种颜色的出现次数
        color_counts = Counter()
        for colors in edge_colors:
            for c in colors:
                color_counts[c] += 1
        info['color_counts'] = dict(color_counts)

        if not color_counts:
            return None, [], info

        # 检查是否交替模式：edges[i] 和 edges[i+2] 应来自同一孩子
        if n >= 4:
            # 提取每条边的 child set（忽略 chromosome）
            edge_children = []
            for colors in edge_colors:
                edge_children.append(set(c[0] for c in colors))
            alternating = all(
                edge_children[i] == edge_children[(i + 2) % n]
                for i in range(n)
            )

            if alternating:
                # 交替块边 → inversion 或 RT/URT
                col_set = set()
                for colors in edge_colors:
                    col_set.update(colors)

                if len(col_set) >= 3:
                    # 3+ 种颜色 → 跨染色体 → RT/URT
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                                      if rare_colors & colors]
                    info['conflict_colors'] = rare_colors
                    etype = 'reciprocal_translocation'
                    return etype, conflict_edges, info
                else:
                    # 2 种颜色 → 染色体内块反转 → inversion
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    if len(rare_colors) > 1 and all(
                            cnt == min_count for cnt in color_counts.values()):
                        rare_colors = {next(iter(color_counts.keys()))}
                    conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                                      if rare_colors & colors]
                    info['conflict_colors'] = rare_colors
                    etype = 'inversion'
                    return etype, conflict_edges, info

        # 3-cycle 在块级可能是 inversion：两条不同孩子颜色的唯一边 + 一条共享边
        if n == 3:
            # 检查：两条唯一边（颜色不同） + 一条共享边
            unique_group = []
            shared_group = []
            for (b1, b2), colors in zip(edges, edge_colors):
                if len(colors) == 1:
                    unique_group.append((b1, b2, next(iter(colors))[0]))
                else:
                    shared_group.append(colors)

            if len(unique_group) == 2 and len(shared_group) == 1:
                # 两条唯一边来自不同孩子 → 块反转
                c1, c2 = unique_group[0][2], unique_group[1][2]
                if c1 != c2:
                    # 3-cycle 中至少 2 个块是多 HOG 块才算 inversion
                    block_sizes = [len(self._blocks.get(bid, [])) for bid in cycle]
                    n_multi = sum(1 for s in block_sizes if s >= 2)
                    if n_multi < 2:
                        # 多为单例块 → 可能为 indel 边界混淆
                        info['pattern'] = 'block_3cycle_inversion_filtered'
                        return None, [], info
                    # 两条唯一边各来自一个孩子 → 2 色交替 = inversion
                    info['pattern'] = 'block_3cycle_inversion'
                    # 根据唯一边确定冲突方：谁是衍生方由外类群判定
                    # 冲突边 = 两条唯一边
                    conflict_edges = [(b1, b2) for b1, b2, _ in unique_group]
                    conflict_colors = {unique_group[0][2]: None, unique_group[1][2]: None}
                    info['conflict_colors'] = conflict_colors
                    etype = 'inversion'
                    return etype, conflict_edges, info
                elif len(unique_group[0]) == 1 and len(unique_group[1]) == 1:
                    # 两条唯一边来自同一个孩子 → 可能是其他事件
                    pass
            elif len(unique_group) == 3:
                # 三条唯一边 → indel 或错误
                etype = 'gene_indel'
                conflict_edges = [(b1, b2) for b1, b2, _ in unique_group]
                info['pattern'] = 'block_3cycle_all_unique'
                return etype, conflict_edges, info

        # 无法分类的块级环
        min_count = min(color_counts.values())
        rare_colors = {c for c, cnt in color_counts.items() if cnt == min_count}
        conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                          if rare_colors & colors]
        info['fallback'] = True
        info['fallback_reason'] = 'unclassified_block_cycle'
        return 'unclassified', conflict_edges, info

    def _resolve_block_structural_events(self, outgroup_adjacency: Set = None):
        """在块级图上检测结构重排事件。

        对每个块级环：
        1. classify_cycle → 判断事件类型
        2. outgroup voting → 确定衍生方
        3. 移除衍生方的块间边
        4. 记录事件

        Returns: 检测到的事件数
        """
        if not hasattr(self, '_block_graph'):
            self._compress_to_block_level()

        iteration = 0
        max_iterations = 100
        events_found = 0

        while iteration < max_iterations:
            iteration += 1
            cycles = self._find_block_cycles()
            if not cycles:
                break
            logger.debug("  [structural] iteration %d: %d cycles found", iteration, len(cycles))

            # Step 1: 对所有环分类
            classified = []
            for cycle in cycles:
                etype, conflict_edges, info = self._classify_block_cycle(cycle)
                logger.debug("  [structural] cycle %s → %s, conflicts=%d, info=%s",
                             [str(b) for b in cycle], etype, len(conflict_edges),
                             {k: v for k, v in info.items() if k != 'edge_colors'})
                if etype is not None and conflict_edges:
                    classified.append((etype, conflict_edges, info, cycle))

            if not classified:
                break

            resolved = 0

            # Step 2: 合并相邻 inversion 3-cycle（同一 inversion 的 N-1 个环）
            inv_cycles = [(ce, info, cyc) for et, ce, info, cyc in classified
                          if et == 'inversion']
            other_cycles = [(et, ce, cyc) for et, ce, info, cyc in classified
                            if et != 'inversion']

            # 合并 inversion: 重叠块的环合并为一个事件
            merged_inv_groups = []
            used_inv = set()
            for i, (ce_i, info_i, cyc_i) in enumerate(inv_cycles):
                if i in used_inv:
                    continue
                group = [i]
                blocks_set = set(cyc_i)
                # 找到所有与当前环共享块的相邻环
                changed = True
                while changed:
                    changed = False
                    for j, (ce_j, info_j, cyc_j) in enumerate(inv_cycles):
                        if j in used_inv or j in group:
                            continue
                        if blocks_set & set(cyc_j):
                            group.append(j)
                            blocks_set.update(cyc_j)
                            changed = True
                used_inv.update(group)
                # 收集组内所有冲突边 + 块
                all_conflict = []
                all_blocks = set()
                for g in group:
                    all_conflict.extend(inv_cycles[g][0])
                    all_blocks.update(inv_cycles[g][2])
                merged_inv_groups.append((all_conflict, all_blocks))

            # Step 3: 处理合并后的 inversion 组 + 其他事件
            for conflict_edges, blocks_set in merged_inv_groups:
                # 去重冲突边
                conflict_edges = list(set(conflict_edges))
                # 外类群投票：合并组只有一个极性
                child_derived = defaultdict(int)
                for (b1, b2) in conflict_edges:
                    colors = self._block_edge_data(b1, b2).get('colors', set())
                    if outgroup_adjacency is not None:
                        hogs1 = self._blocks.get(b1, [])
                        hogs2 = self._blocks.get(b2, [])
                        any_ancestral = False
                        for h1 in hogs1:
                            for h2 in hogs2:
                                if self.has_edge(h1, h2):
                                    u_p = getattr(h1, 'parent',
                                                  h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                    v_p = getattr(h2, 'parent',
                                                  h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                    key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                    if key in outgroup_adjacency:
                                        any_ancestral = True
                                        break
                        for cid, _ in colors:
                            if any_ancestral:
                                pass  # 祖先方不计数
                            else:
                                child_derived[cid] += 1
                    else:
                        for c in colors:
                            child_derived[c[0]] -= 1  # 无外类群回退

                if child_derived:
                    derived_cid = max(child_derived, key=child_derived.get)
                else:
                    # 回退：边上的颜色出现次数较少的
                    cc = Counter()
                    for (b1, b2) in conflict_edges:
                        for c in self._block_edge_data(b1, b2).get('colors', set()):
                            cc[c] += 1
                    derived_cid = min(cc, key=cc.get)[0] if cc else None

                if derived_cid is None:
                    continue

                # 移除衍生方在冲突边上的颜色
                for (b1, b2) in conflict_edges:
                    colors = self._block_edge_data(b1, b2).get('colors', set())
                    remaining = set()
                    for cid, chrom in colors:
                        if cid != derived_cid:
                            remaining.add((cid, chrom))
                    if remaining:
                        self._block_edge_data(b1, b2)['colors'] = remaining
                    else:
                        self._block_remove_edge(b1, b2)

                # 记录合并后的 inversion 事件
                involved_hogs = []
                for bid in blocks_set:
                    involved_hogs.extend(self._blocks.get(bid, []))
                involved_hogs = list(set(involved_hogs))

                event = TAKREvent(
                    event_type='inversion',
                    branch=f"{self.hog_level}-{derived_cid}" if derived_cid else self.hog_level,
                    genes_involved=involved_hogs,
                    desc="inversion (block-level): {} blocks, {} HOGs".format(
                        len(blocks_set), len(involved_hogs)),
                    support=1,
                )
                self.events.append(event)
                resolved += 1
                events_found += 1

                # 同步移除 HOG 级图中衍生方的边
                for (b1, b2) in conflict_edges:
                    for h1 in self._blocks.get(b1, []):
                        for h2 in self._blocks.get(b2, []):
                            if self.has_edge(h1, h2):
                                colors = self.get_colors(h1, h2)
                                for cid, chrom in list(colors):
                                    if cid == derived_cid:
                                        self.remove_edge_color(h1, h2, (cid, chrom))

            # 处理非 inversion 事件（逐环处理）
            for etype, conflict_edges, cycle in other_cycles:
                # gene_indel 和 unclassified 不是真正的结构重排，不移除边
                if etype in ('gene_indel', 'unclassified'):
                    continue
                conflict_edges = list(set(conflict_edges))
                # 检查冲突边是否还在块级图中
                valid_edges = [(b1, b2) for b1, b2 in conflict_edges
                               if self._block_has_edge(b1, b2)]
                if not valid_edges:
                    continue
                # 外类群投票
                child_derived2 = defaultdict(int)
                for (b1, b2) in valid_edges:
                    colors = self._block_edge_data(b1, b2).get('colors', set())
                    if outgroup_adjacency is not None:
                        hogs1 = self._blocks.get(b1, [])
                        hogs2 = self._blocks.get(b2, [])
                        any_ancestral = False
                        for h1 in hogs1:
                            for h2 in hogs2:
                                if self.has_edge(h1, h2):
                                    u_p = getattr(h1, 'parent',
                                                  h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                    v_p = getattr(h2, 'parent',
                                                  h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                    key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                    if key in outgroup_adjacency:
                                        any_ancestral = True
                                        break
                        for cid, _ in colors:
                            if not any_ancestral:
                                child_derived2[cid] += 1
                if child_derived2:
                    derived_cid2 = max(child_derived2, key=child_derived2.get)
                else:
                    cc2 = Counter()
                    for (b1, b2) in valid_edges:
                        for c in self._block_edge_data(b1, b2).get('colors', set()):
                            cc2[c] += 1
                    derived_cid2 = min(cc2, key=cc2.get)[0] if cc2 else None

                if derived_cid2 is None:
                    continue

                # 移除衍生方的边
                for (b1, b2) in valid_edges:
                    colors = self._block_edge_data(b1, b2).get('colors', set())
                    remaining = set()
                    for cid, chrom in colors:
                        if cid != derived_cid2:
                            remaining.add((cid, chrom))
                    if remaining:
                        self._block_edge_data(b1, b2)['colors'] = remaining
                    else:
                        self._block_remove_edge(b1, b2)

                # 记录事件
                involved_hogs2 = []
                for bid in cycle:
                    involved_hogs2.extend(self._blocks.get(bid, []))
                involved_hogs2 = list(set(involved_hogs2))

                event = TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{derived_cid2}" if derived_cid2 else self.hog_level,
                    genes_involved=involved_hogs2,
                    desc="{} (block-level): {} blocks, {} HOGs".format(
                        etype, len(cycle), len(involved_hogs2)),
                    support=1,
                )
                self.events.append(event)
                resolved += 1
                events_found += 1

                # 同步移除 HOG 级边
                for (b1, b2) in valid_edges:
                    for h1 in self._blocks.get(b1, []):
                        for h2 in self._blocks.get(b2, []):
                            if self.has_edge(h1, h2):
                                colors = self.get_colors(h1, h2)
                                for cid, chrom in list(colors):
                                    if cid == derived_cid2:
                                        self.remove_edge_color(h1, h2, (cid, chrom))

            if resolved == 0:
                break

        return events_found

    def _resolve_block_bridge_events(self, outgroup_adjacency: Set = None):
        """在块级图上检测桥接事件。

        块级桥接：块级唯一边连接两个块级共享分量。
        端粒保护：不移除涉及共识端粒块的边。
        """
        if not hasattr(self, '_block_graph'):
            self._compress_to_block_level()

        bg = self._block_graph

        # 块级共享边图
        shared_bg = nx.DiGraph()
        for b1, b2, data in bg.edges(data=True):
            if len(data.get('colors', set())) > 1:
                shared_bg.add_edge(b1, b2)

        if shared_bg.number_of_edges() == 0:
            return 0

        components = list(nx.weakly_connected_components(shared_bg))
        block_to_comp = {}
        for ci, comp in enumerate(components):
            for b in comp:
                block_to_comp[b] = ci

        # 共识端粒块：包含端粒 HOG 的块 (min_children=1 更宽松)
        cons_tels = self.consensus_telomeres(min_children=1)
        telomere_blocks = set()
        for hog in cons_tels:
            bid = self._hog_to_block.get(hog)
            if bid:
                telomere_blocks.add(bid)

        # 共享分量 → 检查是否包含端粒块（真正的染色体）
        # 使用原始共享组件（Phase 2 前），避免 indel 碎片化影响
        if hasattr(self, '_orig_comp_has_telomere'):
            block_to_orig_comp = {}
            for bid, hogs in self._blocks.items():
                for h in hogs:
                    if h in self._hog_to_orig_comp:
                        block_to_orig_comp[bid] = self._hog_to_orig_comp[h]
                        break
            comp_has_telomere = set()
            for bid in block_to_orig_comp:
                if bid in telomere_blocks:
                    comp_has_telomere.add(block_to_orig_comp[bid])
            use_orig_comp = True
        else:
            comp_has_telomere = set()
            for ci, comp in enumerate(components):
                for b in comp:
                    if b in telomere_blocks:
                        comp_has_telomere.add(ci)
                        break
            use_orig_comp = False

        events_found = 0
        n_unique = sum(1 for _, _, d in bg.edges(data=True) if len(d.get('colors', set())) == 1)
        n_shared = sum(1 for _, _, d in bg.edges(data=True) if len(d.get('colors', set())) > 1)
        logger.debug("  [bridge] block graph: %d nodes, %d edges (%d unique, %d shared), %d shared components",
                     bg.number_of_nodes(), bg.number_of_edges(), n_unique, n_shared, len(components))
        logger.debug("  [bridge] telomere blocks: %s, comp_has_telomere: %s",
                     telomere_blocks, comp_has_telomere)
        for b1, b2, data in list(bg.edges(data=True)):
            colors = data.get('colors', set())
            if len(colors) != 1:
                continue
            if use_orig_comp:
                c1 = block_to_orig_comp.get(b1)
                c2 = block_to_orig_comp.get(b2)
            else:
                c1 = block_to_comp.get(b1)
                c2 = block_to_comp.get(b2)
            b1_has_tel = b1 in telomere_blocks
            b2_has_tel = b2 in telomere_blocks
            logger.debug("  [bridge] unique edge %s↔%s: comp=%s↔%s, tel=%s↔%s, same_comp=%s",
                         b1, b2, c1, c2, b1_has_tel, b2_has_tel,
                         c1 == c2 if c1 is not None and c2 is not None else 'N/A')
            if c1 is not None and c2 is not None and c1 != c2:
                # 端粒验证
                has_tel = b1_has_tel or b2_has_tel
                color = next(iter(colors))
                child_id = color[0]

                if outgroup_adjacency is not None:
                    # 检查跨块 HOG 对的 parent 级邻接
                    hogs1 = self._blocks.get(b1, [])
                    hogs2 = self._blocks.get(b2, [])
                    any_ancestral = False
                    matched_key = None
                    for h1 in hogs1:
                        for h2 in hogs2:
                            if self.has_edge(h1, h2):
                                u_p = getattr(h1, 'parent',
                                              h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                v_p = getattr(h2, 'parent',
                                              h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                if key in outgroup_adjacency:
                                    any_ancestral = True
                                    matched_key = key
                                    break
                        if any_ancestral:
                            break
                    logger.debug("  [bridge] %s↔%s: outgroup_check any_ancestral=%s, key=%s",
                                 b1, b2, any_ancestral, matched_key)

                    if any_ancestral:
                        # fission: 外类群有该邻接 → 需要端粒验证防止误报
                        if not has_tel:
                            continue
                        etype = 'fission'
                        suffix = f" (outgroup confirms {child_id} preserved)"
                    else:
                        # 检查 NCF 条件
                        orig_comps = getattr(self, '_original_shared_components', None)
                        if use_orig_comp and orig_comps and c1 < len(orig_comps) and c2 < len(orig_comps):
                            c_size = [len(list(orig_comps[c1])), len(list(orig_comps[c2]))]
                        else:
                            c_size = [len(components[c1]) if c1 < len(components) else 0,
                                      len(components[c2]) if c2 < len(components) else 0]
                        if min(c_size) < 5:
                            etype = 'ncf'
                        else:
                            etype = 'eej'
                        suffix = f" (outgroup: derived in {child_id})"
                else:
                    # 无外类群 (root): 仅在有端粒时记录（保守策略）
                    if not has_tel:
                        continue
                    etype = 'bridge_unclassified'
                    suffix = " (no outgroup, root node)"

                # 确定事件发生的分支:
                # - 唯一边属于 child_id
                # - fission: 另一个孩子丢失了邻接 → 事件在另一个孩子分支
                # - fusion (eej/ncf): child_id 获得了邻接 → 事件在 child_id 分支
                other_children = [c for c in self.children() if c != child_id]
                if etype == 'fission' and other_children:
                    event_branch = f"{self.hog_level}-{other_children[0]}"
                else:
                    event_branch = f"{self.hog_level}-{child_id}"

                # 展开 HOG
                involved_hogs = list(self._blocks.get(b1, [])) + list(self._blocks.get(b2, []))
                involved_hogs = list(set(involved_hogs))

                logger.debug("  [bridge] %s: %s(%d HOGs,comp%d) <-> %s(%d HOGs,comp%d) via %s → branch=%s",
                            etype, b1, len(self._blocks.get(b1, [])), c1,
                            b2, len(self._blocks.get(b2, [])), c2, child_id, event_branch)

                self.events.append(TAKREvent(
                    event_type=etype,
                    branch=event_branch,
                    genes_involved=involved_hogs,
                    desc="{} block-bridge: {} + {}{}".format(
                        etype, b1, b2, suffix),
                    support=1,
                ))

                # 移除边：只移除 fusion (eej/ncf) 的边，保留 fission 边
                # fission: 外类群有该邻接 → 祖先态 → 保留边（另一个孩子丢失了它）
                # fusion: 外类群无该邻接 → 衍生态 → 移除边（这个孩子获得了它）
                if etype in ('eej', 'ncf'):
                    for h1 in self._blocks.get(b1, []):
                        for h2 in self._blocks.get(b2, []):
                            if self.has_edge(h1, h2):
                                self.remove_edge_color(h1, h2, color)
                    bg.remove_edge(b1, b2)
                # fission: 不移除边，只记录事件
                events_found += 1

        return events_found

    # ==================== 桥接冲突检测 (EEJ/NCF) ====================

    def resolve_bridge_events(self, outgroup_adjacency: Set = None):
        """检测桥接冲突: unique edge 连接两个共享边连通分量。

        外类群投票区分 fission vs fusion (EEJ/NCF):
        - 外类群 (parent HOG level) 有 (h1, h2) → 祖先有该连接
          → 保存的孩子有边，丢失的孩子没有 → fission on 丢失方的分支
        - 外类群无 (h1, h2) → 祖先无该连接
          → 有边的孩子是融合方 → EEJ/NCF on 有边方的分支
        - 无外类群 (root): 记录为 bridge_unclassified (需上游处理)

        算法:
        1. 建共享边图 (edges with >1 color → ≥2 children)
        2. 找连通分量
        3. 对每条 unique edge (1 color):
           两端在不同分量 → 桥接冲突
           → 用外类群投票分类 → 移除并记录事件
        """
        # Build shared-edge graph
        shared_G = nx.DiGraph()
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) > 1:  # shared by ≥2 children
                shared_G.add_edge(h1, h2)

        if shared_G.number_of_edges() == 0:
            return

        # Find connected components in shared graph
        components = list(nx.weakly_connected_components(shared_G))
        hog_to_comp = {}
        for ci, comp in enumerate(components):
            for h in comp:
                hog_to_comp[h] = ci

        # Check each unique edge for bridge pattern
        cons_tels = self.consensus_telomeres(min_children=2)

        # 共享分量 → 检查是否包含端粒 HOG（真正的染色体）
        comp_has_telomere = set()
        for ci, comp in enumerate(components):
            for h in comp:
                if h in cons_tels:
                    comp_has_telomere.add(ci)
                    break

        for h1, h2, data in list(self.edges(data=True)):
            if len(data['colors']) != 1:
                continue  # not a unique edge
            c1 = hog_to_comp.get(h1)
            c2 = hog_to_comp.get(h2)
            if c1 is not None and c2 is not None and c1 != c2:
                # 端粒验证：有外类群时不需要，无外类群时要求至少一个分量有端粒
                has_tel = c1 in comp_has_telomere or c2 in comp_has_telomere
                # Bridge: connects two different shared components
                color = next(iter(data['colors']))
                child_id = color[0]
                comp1_size = len(components[c1])
                comp2_size = len(components[c2])

                # === 外类群投票分类 ===
                if outgroup_adjacency is not None:
                    # 将 h1, h2 映射到 parent HOG level 用于外类群查询
                    h1_parent = getattr(h1, 'parent', h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                    h2_parent = getattr(h2, 'parent', h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                    if isinstance(h1_parent, tuple) or isinstance(h1_parent, (int, float)):
                        # HOGrecord 的 parent 是 hog_id string
                        h1_pid = str(h1_parent) if not isinstance(h1_parent, str) else h1_parent
                        h2_pid = str(h2_parent) if not isinstance(h2_parent, str) else h2_parent
                    else:
                        h1_pid = h1_parent
                        h2_pid = h2_parent
                    key = (h1_pid, h2_pid) if h1_pid < h2_pid else (h2_pid, h1_pid)

                    if key in outgroup_adjacency:
                        # fission: 外类群有该邻接 → 需要端粒验证防止误报
                        if not has_tel:
                            continue
                        etype = 'fission'
                        branch_suffix = f" (outgroup confirms {child_id} preserved)"
                    else:
                        # 外类群无该邻接 → 祖先无连接 → 当前孩子融合了它 → EEJ/NCF
                        # NCF: 一端在共享图的端粒位置
                        in_shared_g = h1 in shared_G and h2 in shared_G
                        deg1 = False
                        if in_shared_g:
                            deg1 = shared_G.degree(h1) == 1 or shared_G.degree(h2) == 1
                        if deg1 or min(comp1_size, comp2_size) < 10:
                            etype = 'ncf'
                        else:
                            etype = 'eej'
                        branch_suffix = f" (outgroup: derived in {child_id})"
                else:
                    # 无外类群 (root): 仅在有端粒时记录
                    if not has_tel:
                        continue
                    etype = 'bridge_unclassified'
                    branch_suffix = " (no outgroup, root node)"

                # 确定事件发生的分支
                other_children = [c for c in self.children() if c != child_id]
                if etype == 'fission' and other_children:
                    event_branch = f"{self.hog_level}-{other_children[0]}"
                else:
                    event_branch = f"{self.hog_level}-{child_id}"

                self.events.append(TAKREvent(
                    event_type=etype,
                    branch=event_branch,
                    genes_involved=[h1, h2],
                    desc=f"{etype} bridge: comp{c1} ({comp1_size} HOGs)"
                         f" + comp{c2} ({comp2_size} HOGs)"
                         f" via {child_id} in {color}{branch_suffix}",
                    support=1,
                ))
                self.remove_edge_color(h1, h2, color)

    # ==================== Pipeline postcondition checks ====================

    def _postcondition_linear(self, label: str):
        """Phase 1 postcondition: 每个节点 degree ≤ 2，无并行边。"""
        violations = []
        for n in self.nodes():
            deg = self._degree(n)
            if deg > 2:
                violations.append((n, deg))
        if violations:
            worst = max(violations, key=lambda x: x[1])
            raise RuntimeError(
                f"{label} postcondition failed: {len(violations)} nodes with degree>2, "
                f"worst={worst[0]} deg={worst[1]}")

    def _postcondition_no_cycles(self, label: str):
        """Phase 2/4d postcondition: 图中无环。"""
        cycles = self.find_cycles()
        if cycles:
            raise RuntimeError(
                f"{label} postcondition failed: {len(cycles)} cycles remain, "
                f"first={cycles[0][:5]}...")

    def _postcondition_no_cross_edges(self, label: str):
        """Phase 2 postcondition: 无跨越边（spanning edges）。"""
        spanning = self.find_spanning_edges()
        if spanning:
            raise RuntimeError(
                f"{label} postcondition failed: {len(spanning)} spanning edges remain")

    def _postcondition_all_hogs_assigned(self, label: str):
        """Phase 3 postcondition: 所有 HOG 已分配到块。"""
        if not hasattr(self, '_hog_to_block'):
            raise RuntimeError(f"{label}: _hog_to_block not built")
        unassigned = self.all_hogs() - set(self._hog_to_block.keys())
        if unassigned:
            raise RuntimeError(
                f"{label} postcondition failed: {len(unassigned)} HOGs not in any block")

    def _postcondition_no_bridge_edges(self, label: str):
        """Phase 4f postcondition: 无桥接边。"""
        shared_G = nx.DiGraph()
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) > 1:
                shared_G.add_edge(h1, h2)
        if shared_G.number_of_edges() == 0:
            return
        components = list(nx.weakly_connected_components(shared_G))
        hog_to_comp = {}
        for ci, comp in enumerate(components):
            for h in comp:
                hog_to_comp[h] = ci
        cons_tels = self.consensus_telomeres(min_children=2)
        bridges = []
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) != 1:
                continue
            c1 = hog_to_comp.get(h1)
            c2 = hog_to_comp.get(h2)
            if c1 is not None and c2 is not None and c1 != c2:
                if h1 not in cons_tels and h2 not in cons_tels:
                    bridges.append((h1, h2))
        if bridges:
            raise RuntimeError(
                f"{label} postcondition failed: {len(bridges)} bridge edges remain")

    def _postcondition_all_hogs_covered(self, label: str):
        """Phase 5 postcondition: 所有 HOG 被路径覆盖。"""
        # 需要在 path_cover 之后调用，检查 paths 覆盖所有 HOG
        # 由 resolve_all_events 中直接检查

    # ==================== Phase 1: per-child deduplication ====================

    def _deduplicate_children(self, child_graphs: list, child_source_ids: list,
                              ref_graphs: list = None):
        """Phase 1: 在合图前对每个孩子单独去重。

        用其他孩子+outgroup 作参照，判断重复 HOG 的祖先位置。
        ref_graphs: 参照图列表（其他孩子+outgroup叶图），用于判断祖先态。
        """
        deduped_graphs = []
        for i, (cg, cid) in enumerate(zip(child_graphs, child_source_ids)):
            refs = []
            if ref_graphs:
                refs = [r for j, r in enumerate(ref_graphs) if j != i]
            new_cg = self._deduplicate_single_child(cg, cid, refs)
            deduped_graphs.append(new_cg)
        return deduped_graphs

    def _deduplicate_single_child(self, child_graph, source_id: str,
                                  ref_graphs: list = None):
        """对单个孩子图去重：检测并移除所有重复 HOG。

        不同基因映射到同一 HOG = dup 事件。
        用参照图判断哪个拷贝是祖先态（保留），哪个是衍生态（删除）。

        检测类型：
        - tandem_dup: 连续重复 (pos i 和 i+1 相同 HOG)
        - dispersed_dup: 分散重复 (同染色体或跨染色体)

        用参照图判断祖先拷贝：
        - 某拷贝的邻居在参照图中也是相邻的 → 该拷贝是祖先态
        - 保留祖先态拷贝，删除其他拷贝
        """
        import copy
        new_graph = copy.deepcopy(child_graph)

        ref_graph = self._build_ref_graph(ref_graphs)

        # 优先用 chrom_hogs（保留重复），回退到 chromosomes
        chrom_source = getattr(new_graph, 'chrom_hogs', None)
        if chrom_source:
            chrom_iter = sorted(chrom_source.items())
        else:
            chrom_iter = list(enumerate(new_graph.chromosomes))

        # 收集所有 HOG 出现位置（跨染色体）
        hog_occurrences = defaultdict(list)  # hog_str → [(chrom_idx, pos, hog_obj), ...]
        chrom_hog_lists = {}  # chrom_idx → [hog_obj, ...]
        for chrom_idx, chrom_nodes in chrom_iter:
            hogs = [n for n in chrom_nodes if n not in new_graph.telomeres]
            chrom_hog_lists[chrom_idx] = hogs
            for i, hog in enumerate(hogs):
                hog_occurrences[str(hog)].append((chrom_idx, i, hog))

        # 找所有重复 HOG（出现 >1 次）
        remove_set = set()  # (chrom_idx, pos)
        zero_evidence_hogs = []  # 零证据 HOG 汇总
        for hog_str, occurrences in hog_occurrences.items():
            if len(occurrences) <= 1:
                continue

            # 用参照图最短路径距离判断祖先拷贝
            hog_in_ref = hog_str in ref_graph
            best = None
            best_score = float('inf')    # 越小越好
            all_unreachable = True
            for chrom_idx, pos, hog_obj in occurrences:
                hogs = chrom_hog_lists[chrom_idx]
                score = 0
                reachable = False
                # 左邻居
                if pos > 0:
                    left = str(hogs[pos - 1])
                    if left in ref_graph and hog_str in ref_graph:
                        try:
                            score += nx.shortest_path_length(ref_graph, source=left, target=hog_str)
                            reachable = True
                        except nx.NetworkXNoPath:
                            pass
                # 右邻居
                if pos < len(hogs) - 1:
                    right = str(hogs[pos + 1])
                    if right in ref_graph and hog_str in ref_graph:
                        try:
                            score += nx.shortest_path_length(ref_graph, source=hog_str, target=right)
                            reachable = True
                        except nx.NetworkXNoPath:
                            pass
                if reachable:
                    all_unreachable = False
                if score < best_score:
                    best_score = score
                    best = (chrom_idx, pos, hog_obj)

            if not hog_in_ref or all_unreachable:
                zero_evidence_hogs.append((hog_str, len(occurrences)))

            # 收集非祖先拷贝（基因级）
            for chrom_idx, pos, hog_obj in occurrences:
                if (chrom_idx, pos) == (best[0], best[1]):
                    continue
                remove_set.add((chrom_idx, pos))
                # 暂存基因信息，后续按块分组生成事件
                if not hasattr(self, '_dedup_pending'):
                    self._dedup_pending = []
                self._dedup_pending.append({
                    'source_id': source_id,
                    'chrom_idx': chrom_idx,
                    'pos': pos,
                    'hog_obj': hog_obj,
                    'ancestral_chrom': best[0],
                    'ancestral_pos': best[1],
                    'score': best_score,
                    'hog_str': hog_str,
                })
                logger.debug("  [dedup] %s chrom%d pos%d: %s (ancestral chrom%d pos%d, dist=%d)",
                             source_id, chrom_idx, pos, hog_str,
                             best[0], best[1], best_score)

        if zero_evidence_hogs:
            total_copies = sum(n for _, n in zero_evidence_hogs)
            logger.warning("  [dedup] %s: zero evidence for %d HOGs (%d copies total): %s ...",
                           source_id, len(zero_evidence_hogs), total_copies,
                           ", ".join(h for h, _ in zero_evidence_hogs[:5]))

        # 更新 chrom_hogs：移除被删除的位置
        if remove_set and chrom_source:
            chrom_removals = {}
            for ci, pos in remove_set:
                chrom_removals.setdefault(ci, set()).add(pos)
            for chrom_idx, positions in chrom_removals.items():
                if chrom_idx not in chrom_source:
                    continue
                hogs = chrom_hog_lists.get(chrom_idx, [])
                new_hogs = [h for i, h in enumerate(hogs) if i not in positions]
                old_chrom = chrom_source[chrom_idx]
                new_chrom = []
                hog_idx = 0
                for node in old_chrom:
                    if node in new_graph.telomeres:
                        new_chrom.append(node)
                    elif hog_idx < len(new_hogs):
                        new_chrom.append(new_hogs[hog_idx])
                        hog_idx += 1
                chrom_source[chrom_idx] = new_chrom

        # 按块分组生成事件（连续位置合并为一个事件）
        self._flush_dedup_events(source_id)

        return new_graph

    @staticmethod
    def _build_ref_graph(ref_graphs):
        """从参照图列表构建无向图（最短路径距离查询）。"""
        g = nx.Graph()
        if not ref_graphs:
            return g
        for ref in ref_graphs:
            ch_source = getattr(ref, 'chrom_hogs', None)
            if ch_source:
                for ci, hogs in ch_source.items():
                    gene_hogs = [n for n in hogs if n not in ref.telomeres]
                    for i in range(len(gene_hogs) - 1):
                        a, b = str(gene_hogs[i]), str(gene_hogs[i + 1])
                        g.add_edge(a, b)
            else:
                for chrom in ref.chromosomes:
                    gene_hogs = [n for n in chrom if n not in ref.telomeres]
                    for i in range(len(gene_hogs) - 1):
                        a, b = str(gene_hogs[i]), str(gene_hogs[i + 1])
                        g.add_edge(a, b)
        return g

    def _flush_dedup_events(self, source_id):
        if not hasattr(self, '_dedup_pending') or not self._dedup_pending:
            return
        from itertools import groupby
        dedup_gene_count = 0
        tandem_blocks = 0
        # 按 (source_id, chrom_idx) 分组
        pending = sorted(self._dedup_pending,
                         key=lambda x: (x['source_id'], x['chrom_idx'], x['pos']))
        for (sid, ci), group in groupby(
                pending, key=lambda x: (x['source_id'], x['chrom_idx'])):
            items = sorted(group, key=lambda x: x['pos'])
            # 连续位置合并为块
            blocks = []
            current_block = [items[0]]
            for item in items[1:]:
                if item['pos'] == current_block[-1]['pos'] + 1:
                    current_block.append(item)
                else:
                    blocks.append(current_block)
                    current_block = [item]
            blocks.append(current_block)

            for block in blocks:
                genes = [it['hog_obj'] for it in block]
                if len(block) == 1:
                    extra = len(block)
                    anc_chrom = block[0].get('ancestral_chrom')
                    if anc_chrom is not None and ci == anc_chrom:
                        event_type = 'dispersed_dup_intra'
                    else:
                        event_type = 'dispersed_dup_inter'
                else:
                    event_type = 'tandem_dup'
                    extra = len(block) - 1
                    tandem_blocks += 1
                dedup_gene_count += extra
                self.events.append(TAKREvent(
                    event_type=event_type,
                    branch=f"{self.hog_level}-{sid}",
                    genes_involved=genes,
                    desc=f"{event_type}: {extra} extra copies at chrom{ci} "
                         f"pos{block[0]['pos']}-{block[-1]['pos']}",
                    support=extra,
                ))
                logger.debug("  [dedup] %s: %s block of %d genes (%d extra) at chrom%d pos%d-%d",
                            sid, event_type, len(genes), extra, ci,
                            block[0]['pos'], block[-1]['pos'])
        if dedup_gene_count + tandem_blocks != len(self._dedup_pending):
            logger.warning("  [dedup] %s: removed %d positions != extra %d + tandem %d = %d",
                           source_id, len(self._dedup_pending),
                           dedup_gene_count, tandem_blocks,
                           dedup_gene_count + tandem_blocks)
        else:
            logger.info("  [dedup] %s: %d = %d + %d ✓",
                        source_id, len(self._dedup_pending),
                        dedup_gene_count, tandem_blocks)
        self._dedup_pending = []

    def resolve_all_events(self, outgroup_graph = None,
                           hog_parent = None,
                           min_hogs: int = 3,
                           gfa_debug: bool = False,
                           gfa_prefix: str = None) -> List[TAKREvent]:
        """Telomere-centric event resolution pipeline.

        Pipeline (简单→复杂，每步检查 postcondition，失败则报错停止):
        ─────────────────────────────────────────────────────────────
        1. 每个孩子内部: duplication resolve (tandem, dispersed, proximal, seg_dup)
           Postcondition: 每个孩子图线性 (deg<=2, no parallel edges)
        2. 合图 + 单基因 indel/loss/gain (HOG 级)
           Postcondition: 无跨越边, 无环, 无分支
        3. 共线性块压缩
           Postcondition: 所有 HOG 已分配到块
        4a. deletion / insertion (单基因+多基因 indel)
           Postcondition: 无大段不对称
        4b. inversion (方向冲突)
        4c. unidir_trans (单向转移)
        4d. internal_inversion (内部倒位)
        4e. RT / URT (相互易位)
        4f. EEJ / NCF / fission (桥接事件，最后)
        5. 端粒约束路径覆盖
           Postcondition: 所有 HOG 覆盖, 每条染色体 2 个端粒
        ─────────────────────────────────────────────────────────────
        """
        def _gfa_out(phase_label, use_blocks=True):
            """输出 GFA 调试文件。"""
            if not gfa_debug or not gfa_prefix:
                return
            suffix = ".hog.gfa" if not use_blocks else ".gfa"
            path = f"{gfa_prefix}.{phase_label}{suffix}"
            with open(path, 'w') as fout:
                fout.write(f"H\tphase:{phase_label}\tlevel:{self.hog_level}\t"
                           f"nodes:{self.node_count()}\tedges:{self.edge_count()}\n")
                self.to_gfa(fout, use_blocks=use_blocks)
            logger.info("  [gfa] %s -> %s", phase_label, path)

        n_events_before = len(self.events)

        # Legacy: existing sub-methods use outgroup_adjacency/outgroups dicts.
        # Not updated yet — set to None (they fall back to no-outgroup).
        outgroup_adjacency = None
        outgroups = None

        # Phase 1 dedup 在 orchestrator 中 add_child 前完成

        _gfa_out("p1_merged", use_blocks=False)  # HOG 级
        _gfa_out("p1_merged")

        # 保存原始共享组件（用于 bridge 检测验证）
        self._save_original_shared_components()

        # ====== Phase 2: 方向调和 ======
        self.harmonize_directions()
        logger.info("  [colored] Phase 2 (harmonize): done, %d nodes, %d edges",
                     self.node_count(), self.edge_count())
        _gfa_out("p2_harmonize")

        # ====== Phase 3: 共线性块压缩 ======
        self._build_synteny_blocks()
        self._compress_to_block_level()
        logger.info("  [colored] Phase 3 (blocks): %d blocks (%.1f HOG/block avg)",
                     len(self._blocks),
                     self.node_count() / max(len(self._blocks), 1))
        try:
            self._postcondition_all_hogs_assigned("Phase 3")
        except RuntimeError as e:
            logger.warning("  [colored] %s", e)

        _gfa_out("p3_blocks")

        # ====== Phase 4: 块级事件 (按复杂度递增) ======

        # Phase 4a: deletion / insertion (单基因+多基因 indel)
        n_seg = self._resolve_seg_events(outgroup_graph=outgroup_graph,
                                           hog_parent=hog_parent)
        if n_seg:
            from collections import defaultdict
            for etype in ('deletion', 'insertion', 'reciprocal_indel'):
                seg_events = [e for e in self.events if e.event_type == etype]
                if not seg_events:
                    continue
                child_counts = defaultdict(int)
                child_lens = defaultdict(list)
                for e in seg_events:
                    cid = e.branch.split('-')[-1] if '-' in e.branch else e.branch
                    child_counts[cid] += 1
                    child_lens[cid].append(e.support)
                for cid in sorted(child_counts):
                    lens = child_lens[cid]
                    len_str = f"len {min(lens)}-{max(lens)}" if len(lens) > 1 else f"len {lens[0]}"
                    logger.info("  [Phase 4a] %s events: %s=%d (%s)",
                                cid, etype, child_counts[cid], len_str)
        else:
            logger.warning("  [colored] Phase 4a: 0 deletion events — unexpected")

        _gfa_out("p4a_seg")

        # Phase 4b: inversion (方向冲突 = 真倒位，调和后)
        n_inv_before = len(self.events)
        n_inv = self._detect_inversions()
        n_inv_events = len(self.events) - n_inv_before
        if 1:
            inv_types = Counter(e.event_type for e in self.events[n_inv_before:])
            inv_str = ", ".join(f"{t}={c}" for t, c in sorted(inv_types.items()))
            logger.info("  [colored] Phase 4b (inversions): %d events [%s]", n_inv_events, inv_str)

        _gfa_out("p4b_inv")

        # Phase 4c: unidir_trans
        n_ut = self._resolve_unidir_trans(outgroup_adjacency=outgroup_adjacency)
        if 1:
            ut_events = [e for e in self.events if e.event_type == 'unidir_trans']
            ut_str = ", ".join(f"{e.branch}" for e in ut_events[-5:])
            logger.info("  [colored] Phase 4b (unidir_trans): %d events [%s]", n_ut, ut_str)

        _gfa_out("p4c_ut")

        # ====== Phase 4c-4e: 块级结构重排 (inversion, RT) ======
        n_before_struct = len(self.events)
        n2, e2 = self.node_count(), self.edge_count()
        n_struct = self._resolve_block_structural_events(
            outgroup_adjacency=outgroup_adjacency)
        n3, e3 = self.node_count(), self.edge_count()
        n_struct_events = len(self.events) - n_before_struct
        struct_types = Counter(e.event_type for e in self.events[n_before_struct:])
        struct_str = ", ".join(f"{t}={c}" for t, c in sorted(struct_types.items())) if struct_types else "none"
        logger.info("  [colored] Phase 4c-4e (structural): %d events [%s], nodes %d→%d (-%d), edges %d→%d (-%d)",
                     n_struct_events, struct_str,
                     n2, n3, n2 - n3, e2, e3, e2 - e3)

        _gfa_out("p4e_struct")

        # ====== Phase 4f: 块级桥接 (EEJ/NCF/fission) ======
        n_before_bridge = len(self.events)
        n4, e4 = self.node_count(), self.edge_count()
        n_bridge = self._resolve_block_bridge_events(
            outgroup_adjacency=outgroup_adjacency)
        n5, e5 = self.node_count(), self.edge_count()
        n_bridge_events = len(self.events) - n_before_bridge
        bridge_types = Counter(e.event_type for e in self.events[n_before_bridge:])
        bridge_str = ", ".join(f"{t}={c}" for t, c in sorted(bridge_types.items())) if bridge_types else "none"
        logger.info("  [colored] Phase 4f (bridge): %d events [%s], nodes %d→%d (-%d), edges %d→%d (-%d)",
                     n_bridge_events, bridge_str,
                     n4, n5, n4 - n5, e4, e5, e4 - e5)
        try:
            self._postcondition_no_bridge_edges("Phase 4f")
        except RuntimeError as e:
            logger.warning("  [colored] %s", e)

        _gfa_out("p4f_bridge")

        # ====== Phase 5: 端粒约束路径覆盖 ======
        paths = self.path_cover()
        self._cached_paths = paths  # 缓存给 to_ancestral_graph 复用
        all_hogs = self.all_hogs()
        covered = set()
        for p in paths:
            covered.update(p)
        uncovered = all_hogs - covered
        if uncovered:
            logger.warning("  [colored] Phase 5: %d HOGs not covered by path_cover",
                           len(uncovered))
        n_chrom = len(paths)
        logger.info("  [colored] Phase 5 (path cover): %d chromosomes", n_chrom)

        _gfa_out("p5_paths")

        # 过滤小事件
        kept = []
        for e in self.events:
            if e.event_type in ('eej', 'ncf', 'fission', 'bridge_unclassified',
                                'reciprocal_translocation',
                                'unbalanced_reciprocal_translocation',
                                'inversion', 'telomere_inversion'):
                kept.append(e)
            elif len(e.genes_involved) >= min_hogs:
                kept.append(e)
        self.events = kept

        type_counts = Counter(e.event_type for e in self.events)
        type_str = ', '.join(f"{t}={c}" for t, c in sorted(type_counts.items()))
        logger.info("  [colored] done: %d events %s (after min_hogs=%d), %d chroms",
                     len(self.events), type_str, min_hogs, n_chrom)

        # === 染色体数一致性检查 (按孩子分支) ===
        # 每个孩子独立推算祖先染色体数:
        #   祖先_chroms = child_chroms - fissions + fusions
        # 所有孩子推算结果应一致
        new_events = self.events[n_events_before:]
        # 按分支分组事件
        branch_events = defaultdict(lambda: {'fission': 0, 'fusion': 0})
        for e in new_events:
            # branch 格式: "N1-Sp_1"
            branch = e.branch
            if e.event_type == 'fission':
                branch_events[branch]['fission'] += 1
            elif e.event_type in ('eej', 'ncf'):
                branch_events[branch]['fusion'] += 1

        inferred_ancestor_chroms = {}
        for cid, child_chroms in self._child_chrom_counts.items():
            branch_key = f"{self.hog_level}-{cid}"
            ev = branch_events.get(branch_key, {'fission': 0, 'fusion': 0})
            # 祖先 = 孩子 - fissions + fusions
            anc_chroms = child_chroms - ev['fission'] + ev['fusion']
            inferred_ancestor_chroms[cid] = anc_chroms
            logger.info("  [chrom] %s: %d chroms, fission=%d, fusion=%d → ancestor=%d",
                        cid, child_chroms, ev['fission'], ev['fusion'], anc_chroms)

        # 检查所有孩子推算的祖先染色体数是否一致
        values = list(inferred_ancestor_chroms.values())
        if not values:
            pass  # no children, skip consistency check
        elif len(set(values)) > 1:
            logger.error("  [chrom] INCONSISTENT ancestor chroms: %s", inferred_ancestor_chroms)
        else:
            inferred = values[0]
            logger.info("  [chrom] consistent: ancestor inferred %d chromosomes", inferred)
            # 对比 path_cover 实际结果
            if inferred != n_chrom:
                logger.error("  [chrom] MISMATCH: inferred=%d, actual path_cover=%d",
                             inferred, n_chrom)
            else:
                logger.info("  [chrom] OK: inferred=%d == path_cover=%d",
                            inferred, n_chrom)

        return self.events

    # ==================== Phase 4a/4b: seg events & unidir_trans ====================

    def _save_original_shared_components(self):
        """保存 Phase 2 前的原始共享组件和端粒信息。"""
        import networkx as nx
        shared_G = nx.DiGraph()
        for h1, h2, data in self.edges(data=True):
            if len(data['colors']) >= 2:
                shared_G.add_edge(h1, h2)

        self._original_shared_components = list(nx.weakly_connected_components(shared_G))
        hog_to_orig_comp = {}
        for ci, comp in enumerate(self._original_shared_components):
            for h in comp:
                hog_to_orig_comp[h] = ci
        self._hog_to_orig_comp = hog_to_orig_comp

        cons_tels = self.consensus_telomeres(min_children=1)
        orig_comp_has_telomere = set()
        for ci, comp in enumerate(self._original_shared_components):
            for h in comp:
                if h in cons_tels:
                    orig_comp_has_telomere.add(ci)
                    break
        self._orig_comp_has_telomere = orig_comp_has_telomere

        logger.debug("  [orig] %d shared components, %d with telomere",
                     len(self._original_shared_components),
                     len(orig_comp_has_telomere))

    def _score_block_edge(self, outgroup_graph, ga, gb, hog_parent=None):
        """ga 尾部 HOG 与 gb 首部 HOG 在外群图中是否存在边。
        
        通过 hog_parent 将 HOG 映射到外群图同层级再查边。
        """
        if outgroup_graph is None or not hasattr(self, '_blocks'):
            return 0
        ha_list = self._blocks.get(ga, [])
        hb_list = self._blocks.get(gb, [])
        if not ha_list or not hb_list:
            return 0
        ha = str(ha_list[-1])
        hb = str(hb_list[0])
        if hog_parent:
            ha = hog_parent.get(ha, ha)
            hb = hog_parent.get(hb, hb)
        return 1 if outgroup_graph.has_edge(ha, hb) else 0

    def _score_block_path(self, outgroup_graph, block_ids, hog_parent=None):
        """路径中相邻块之间在外群图中存在的边数。遍历 zip 逐对调用 _score_block_edge。"""
        if outgroup_graph is None or len(block_ids) < 2:
            return 0
        return sum(self._score_block_edge(outgroup_graph, a, b, hog_parent)
                   for a, b in zip(block_ids, block_ids[1:]))

    def _resolve_seg_events(self, outgroup_graph=None, hog_parent=None) -> int:
        """Phase 4a: 气泡结构驱动的 deletion / insertion / reciprocal_indel。

        按 (N1, N2) 邻居对分组单物种块：
        - 0 块：跳过
        - 1 块 B（child X）：简单 indel，有 shortcut 边才处理
        - 2 块 Bx、By（X≠Y）：reciprocal_indel
        外群边打分定祖先态，高分路径保留，低分删除。
        """
        if not hasattr(self, '_blocks') or not self._blocks:
            return 0

        bg = self._block_graph
        if not bg:
            return 0

        n_blocks_removed = 0
        n_edges_removed = 0
        n_insertions = 0
        n_deletions = 0
        n_fallback = 0
        n_simple = 0
        n_reciprocal = 0
        n_reciprocal_resolved = 0
        removed_blocks = set()

        pre_nodes = bg.number_of_nodes()
        pre_edges = bg.number_of_edges()

        n_single_species = sum(
            1 for bid in self._blocks
            if len(set(c for c, _ in bg.nodes[bid].get('sources', set()))) == 1
        )

        # === 按 (N1, N2) 邻居对分组单物种块 ===
        from collections import defaultdict
        groups = defaultdict(list)  # (n1, n2) → [(bid, cid, hogs)]

        for bid in list(self._blocks.keys()):
            if bid in removed_blocks:
                continue
            sources = bg.nodes[bid].get('sources', set())
            owner_set = set(c for c, _ in sources)
            if len(owner_set) != 1:
                continue
            own_cid = next(iter(owner_set))
            preds = [p for p in bg.predecessors(bid)
                     if any(c == own_cid for c, _ in
                            self._block_edge_data(p, bid).get('colors', set()))]
            succs = [s for s in bg.successors(bid)
                     if any(c == own_cid for c, _ in
                            self._block_edge_data(bid, s).get('colors', set()))]
            if len(preds) != 1 or len(succs) != 1:
                continue
            n1, n2 = preds[0], succs[0]
            n1_srcs = set(c for c, _ in bg.nodes[n1].get('sources', set()))
            n2_srcs = set(c for c, _ in bg.nodes[n2].get('sources', set()))
            if len(n1_srcs) < 2 or len(n2_srcs) < 2:
                continue
            groups[(n1, n2)].append((bid, own_cid, self._blocks[bid]))

        n_groups = len(groups)
        logger.debug("  [bubble] %d single-species blocks → %d neighbor-pair groups",
                     sum(len(v) for v in groups.values()), n_groups)

        for (n1, n2), items in groups.items():
            """
            对每对共享邻居 (n1, n2)，中间有 1~2 条候选路径（气泡结构）：
            
            简单 indel（1 个单物种块 B）:
              path1 = N1→B→N2  (own_cid 的块路径)
              path2 = N1→N2    (other_cid 的 shortcut 边，如果存在)
            
            reciprocal_indel（2 个单物种块）:
              path1 = N1→B1→N2 (c1 的块路径)
              path2 = N1→B2→N2 (c2 的块路径)

            处理流程:
              1. 显式构造 path1, path2（含块/边、孩子、HOG 列表）
              2. 打分：_score_block_path / _score_block_edge
              3. 选 winner：高分胜出；同分取块路径优先、更长的优先
              4. 加颜色：winner 边加上 loser 的孩子颜色（统一为共享边）
              5. 删 loser：块→remove_node；shortcut 边→remove_edge
            """

            # ── Step 1: 构造 path1 ──
            bid, cid, hogs = items[0]
            path1_bid = bid
            path1_cid = cid
            path1_hogs = hogs
            path1_edges = [(n1, bid), (bid, n2)]

            # ── Step 2: 构造 path2 ──
            # 情况 A: reciprocal_indel，第二个块
            # 情况 B: 简单 indel，shortcut 边作为 path2
            path2_bid = None
            path2_cid = None
            path2_hogs = []
            path2_edges = []

            if len(items) == 2:
                # reciprocal: 两个孩子各有一个平行块
                p2_bid, p2_cid, p2_hogs = items[1]
                if p2_cid == cid:
                    continue  # 同一孩子，跳过
                path2_bid = p2_bid
                path2_cid = p2_cid
                path2_hogs = p2_hogs
                path2_edges = [(n1, p2_bid), (p2_bid, n2)]
                is_simple = False
            elif bg.has_edge(n1, n2):
                # 简单 indel: shortcut 边作为 path2
                other_cid = None
                for c in self.children():
                    if c != cid:
                        other_cid = c
                        break
                if other_cid is None:
                    continue
                # 确认 shortcut 带 other_cid 颜色
                sc_colors = self._block_edge_data(n1, n2).get('colors', set())
                if not any(c == other_cid for c, _ in sc_colors):
                    continue
                path2_cid = other_cid
                path2_edges = [(n1, n2)]
                is_simple = True
            else:
                continue

            # ── Step 3: 打分 ──
            # 块路径用 _score_block_path(块序列)，shortcut 用 _score_block_edge(单边)
            score1 = self._score_block_path(outgroup_graph,
                                            [n1, path1_bid, n2], hog_parent)
            if path2_bid is not None:
                score2 = self._score_block_path(outgroup_graph,
                                                [n1, path2_bid, n2], hog_parent)
            else:
                score2 = self._score_block_edge(outgroup_graph, n1, n2, hog_parent)

            # ── Step 4: 选 winner ──
            has_og = outgroup_graph is not None
            if score1 > score2:
                winner_bid, winner_cid, winner_hogs, winner_edges = \
                    path1_bid, path1_cid, path1_hogs, path1_edges
                loser_bid, loser_cid, loser_hogs = \
                    path2_bid, path2_cid, path2_hogs
            elif score2 > score1:
                winner_bid, winner_cid, winner_hogs, winner_edges = \
                    path2_bid, path2_cid, path2_hogs, path2_edges
                loser_bid, loser_cid, loser_hogs = \
                    path1_bid, path1_cid, path1_hogs
            else:
                # 同分：块路径优先；都块→更长的优先；都 shortcut→path1
                if path1_bid is not None:
                    if path2_bid is not None:
                        # 两块路径，取更长的
                        if len(path1_hogs) >= len(path2_hogs):
                            winner_bid, winner_cid, winner_hogs, winner_edges = \
                                path1_bid, path1_cid, path1_hogs, path1_edges
                            loser_bid, loser_cid, loser_hogs = \
                                path2_bid, path2_cid, path2_hogs
                        else:
                            winner_bid, winner_cid, winner_hogs, winner_edges = \
                                path2_bid, path2_cid, path2_hogs, path2_edges
                            loser_bid, loser_cid, loser_hogs = \
                                path1_bid, path1_cid, path1_hogs
                    else:
                        # path1 是块，path2 是 shortcut → 块优先
                        winner_bid, winner_cid, winner_hogs, winner_edges = \
                            path1_bid, path1_cid, path1_hogs, path1_edges
                        loser_bid, loser_cid, loser_hogs = \
                            path2_bid, path2_cid, path2_hogs
                elif path2_bid is not None:
                    # path2 是块，path1 是 shortcut → 块优先
                    winner_bid, winner_cid, winner_hogs, winner_edges = \
                        path2_bid, path2_cid, path2_hogs, path2_edges
                    loser_bid, loser_cid, loser_hogs = \
                        path1_bid, path1_cid, path1_hogs
                else:
                    # 两条都是 shortcut → path1 胜
                    winner_bid, winner_cid, winner_hogs, winner_edges = \
                        path1_bid, path1_cid, path1_hogs, path1_edges
                    loser_bid, loser_cid, loser_hogs = \
                        path2_bid, path2_cid, path2_hogs

            # ── Step 5: 加颜色 ──
            # winner 路径保留为祖先态，补上 loser 的孩子颜色
            # 使祖先路径成为共享边（两边孩子都有）
            for u, v in winner_edges:
                wdata = self._block_edge_data(u, v)
                wcols = wdata.get('colors', set())
                ch = next((s[1] for s in bg.nodes[n1].get('sources', set())
                           if isinstance(s, tuple) and len(s) == 2
                           and s[0] == loser_cid), 0)
                wcols.add((loser_cid, ch))
                wdata['colors'] = wcols

            # ── Step 6: 删 loser ──
            if loser_bid is not None:
                # loser 是块 → 删块节点和所有邻接边
                # （祖先态不需要这个块，它是衍生的插入）
                n_edges_removed += bg.degree(loser_bid)
                bg.remove_node(loser_bid)
                del self._blocks[loser_bid]
                removed_blocks.add(loser_bid)
                n_blocks_removed += 1
                n_insertions += 1
                etype = 'insertion' if is_simple else 'reciprocal_indel'
            else:
                # loser 是 shortcut 边 → 删边
                # （衍生态 shortcut，祖先态走 winner 块路径）
                if bg.has_edge(n1, n2):
                    bg.remove_edge(n1, n2)
                    n_edges_removed += 1
                n_deletions += 1
                etype = 'deletion'

            self.events.append(TAKREvent(
                event_type=etype,
                branch=f"{self.hog_level}-{loser_cid}",
                genes_involved=list(loser_hogs) if loser_hogs else [],
                desc=f"{etype}: {'block ' + str(loser_bid) if loser_bid else 'shortcut ' + str((n1, n2))}"
                     f" resolved ({n1}↔{n2}, winner={winner_cid}, "
                     f"scores={score1}/{score2})",
                support=len(loser_hogs) if loser_hogs else 1,
            ))
            if is_simple:
                n_simple += 1
            else:
                n_reciprocal += 1
                n_reciprocal_resolved += 1

        # 统计 loser 反推 winner 得 score 分布
        # （保留原 simple_dist / score_dist 逻辑，从 event 信息反填）

        if n_blocks_removed or n_edges_removed:
            post_nodes = bg.number_of_nodes()
            post_edges = bg.number_of_edges()
            logger.info("  [seg] removed: %d blocks, %d edges "
                        "(%d insertions, %d deletions%s)",
                        n_blocks_removed, n_edges_removed,
                        n_insertions, n_deletions,
                        f", {n_fallback} fallback" if n_fallback else "")
            logger.info("  [seg] graph: %d→%d nodes, %d→%d edges (Δn=%+d, Δe=%+d)",
                        pre_nodes, post_nodes, pre_edges, post_edges,
                        post_nodes - pre_nodes, post_edges - pre_edges)
        logger.info("  [seg] single-species blocks: %d total → %d resolved as indel"
                    " (%d simple + %d reciprocal%s)",
                    n_single_species,
                    n_simple + n_reciprocal_resolved,
                    n_simple, n_reciprocal,
                    f", {n_reciprocal_resolved}/{n_reciprocal} resolved"
                    if n_reciprocal else "")


        return n_blocks_removed + n_deletions

    def _resolve_unidir_trans(self, outgroup_adjacency: Set = None) -> int:
        """Phase 4b: 检测 unidirectional translocation。

        模式: A 有 a-b-c-d，B 有 a-x-d (b,c 被 x 替换，x 只在 B 中)。
        在块级图中: 两个孩子共享 a 和 d，但中间路径不同。
        """
        # 简化实现: 检查块级图中的非循环结构差异
        # 完整实现需要对每个孩子的块序列做比对
        # 当前返回 0，后续补充
        return 0

    # ==================== 可视化 ====================

    def draw_block_graph(self, outpath: str, dpi: int = 200,
                         title: str = None, show_hogs: bool = False):
        """Draw block-level colored adjacency graph using Graphviz.

        Nodes = synteny blocks (sized by HOG count).
        Edges:
          - Shared (multi-child): thick black solid
          - Unique (single-child): thin colored solid, one color per child
        Telomere-adjacent blocks: double-octagon shape, orange fill.

        Args:
            outpath: Output PNG path
            dpi: Resolution
            title: Plot title
            show_hogs: If True, show individual HOG names in block labels
        """
        try:
            import graphviz
            # Check if dot binary is available
            import shutil
            if not shutil.which('dot'):
                raise ImportError('dot binary not found')
        except (ImportError, Exception):
            logger.debug("graphviz dot not available, using matplotlib fallback")
            self._draw_block_graph_mpl(outpath, dpi, title)
            return

        if not hasattr(self, '_block_graph') or not self._block_graph:
            self._build_synteny_blocks()
            self._compress_to_block_level()

        bg = self._block_graph
        if bg.number_of_nodes() == 0:
            return

        # Build telomere block set
        telomere_blocks = set()
        for cid, tels in self._child_telomeres.items():
            for tel_hog in tels:
                bid = self._hog_to_block.get(tel_hog)
                if bid:
                    telomere_blocks.add(bid)

        # Cycle nodes
        cycles = []
        try:
            cycles = nx.cycle_basis(bg)
        except Exception:
            pass
        cycle_nodes = set()
        for cyc in cycles:
            cycle_nodes.update(cyc)

        # Child color palette (saturated, distinct)
        children = sorted(self.children())
        palette = ['#c0392b', '#2980b9', '#27ae60', '#f39c12',
                   '#8e44ad', '#16a085', '#d35400', '#2c3e50']
        child_colors = {c: palette[i % len(palette)]
                        for i, c in enumerate(children)}

        # Build graphviz DOT
        dot = graphviz.Digraph(format='png')
        dot.attr(rankdir='LR', label=title or f'Block Graph: {self.hog_level}',
                 labelloc='t', fontsize='18', fontname='Helvetica-Bold',
                 fontcolor='#1a1a2e', bgcolor='white', pad='0.8',
                 nodesep='0.8', ranksep='1.2', ratio='compress', size='16,12',
                 splines='polyline')
        dot.attr('node', fontname='Helvetica-Bold', fontsize='11',
                 fontcolor='#1a1a2e', margin='0.15,0.08')
        dot.attr('edge', fontname='Helvetica', fontsize='8')

        # Nodes
        for n in bg.nodes():
            hog_list = self._blocks.get(n, [])
            hog_count = len(hog_list)
            label = f'{n}\n({hog_count})'
            if show_hogs and hog_count <= 5:
                label += '\n' + ', '.join(str(h) for h in hog_list[:5])

            attrs = {'label': label,
                     'width': str(round(max(0.6, math.log2(hog_count + 1) * 0.4), 2)),
                     'height': str(round(max(0.4, math.log2(hog_count + 1) * 0.3), 2)),
                     'fontcolor': '#1a1a2e', 'fontsize': '14'}

            if n in telomere_blocks:
                attrs['shape'] = 'doubleoctagon'
                attrs['style'] = 'filled,bold'
                attrs['fillcolor'] = '#fff3cd'  # light yellow
                attrs['color'] = '#e67e22'       # orange border
                attrs['penwidth'] = '2.5'
                attrs['fontcolor'] = '#7d5a00'   # dark gold text
            elif n in cycle_nodes:
                attrs['shape'] = 'octagon'
                attrs['style'] = 'filled,bold'
                attrs['fillcolor'] = '#fce4ec'   # light pink
                attrs['color'] = '#c62828'       # dark red border
                attrs['penwidth'] = '2.5'
                attrs['fontcolor'] = '#b71c1c'   # dark red text
            else:
                attrs['shape'] = 'box'
                attrs['style'] = 'filled,rounded'
                attrs['fillcolor'] = '#e3f2fd'   # light blue
                attrs['color'] = '#1565c0'       # blue border
                attrs['penwidth'] = '1.5'

            dot.node(str(n), **attrs)

        # Edges
        for h1, h2, data in bg.edges(data=True):
            colors = data.get('colors', set())
            child_ids = sorted(set(c for c, _ in colors))

            if len(child_ids) > 1:
                # Shared edge: dark navy, very thick
                dot.edge(str(h1), str(h2),
                         color='#1a1a2e', penwidth='3.5',
                         style='solid')
            elif len(child_ids) == 1:
                cid = child_ids[0]
                color = child_colors.get(cid, '#7f8c8d')
                dot.edge(str(h1), str(h2),
                         color=color, penwidth='2.0',
                         style='solid')

        # Legend as subgraph
        with dot.subgraph(name='cluster_legend') as legend:
            legend.attr(label='Legend', style='dashed', fontsize='12',
                        fontcolor='#1a1a2e', bgcolor='#f5f5f5')
            legend.node('_lg_shared', 'shared\n(multi-child)',
                        shape='box', style='filled,rounded',
                        fillcolor='#e3f2fd', color='#1565c0',
                        fontcolor='#1a1a2e')
            legend.node('_lg_cycle', 'in cycle',
                        shape='octagon', style='filled,bold',
                        fillcolor='#fce4ec', color='#c62828',
                        fontcolor='#b71c1c')
            legend.node('_lg_tel', 'telomere',
                        shape='doubleoctagon', style='filled,bold',
                        fillcolor='#fff3cd', color='#e67e22',
                        fontcolor='#7d5a00')
            for cid in children:
                color = child_colors.get(cid, '#7f8c8d')
                legend.node(f'_lg_{cid}', f'{cid}\n(unique)',
                            shape='box', style='filled,rounded',
                            fillcolor=color, fontcolor='white')
            legend.edge('_lg_shared', '_lg_cycle', style='invis')

        # Render
        base = outpath.rsplit('.', 1)[0]
        dot.render(base, cleanup=True)
        logger.info("  [viz] block graph saved to %s", base + '.png')

    def _draw_block_graph_mpl(self, outpath: str, dpi: int = 200,
                              title: str = None):
        """Fallback matplotlib block graph (original implementation)."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not available, skipping block graph")
            return

        if not hasattr(self, '_block_graph') or not self._block_graph:
            self._build_synteny_blocks()
            self._compress_to_block_level()

        bg = self._block_graph
        if bg.number_of_nodes() == 0:
            return

        fig, ax = plt.subplots(1, 1, figsize=(max(16, bg.number_of_nodes() * 0.5), 10))
        try:
            pos = nx.kamada_kawai_layout(bg)
        except Exception:
            pos = nx.spring_layout(bg, k=3.0, seed=42, iterations=100)

        node_sizes = []
        for n in bg.nodes():
            hog_count = len(self._blocks.get(n, []))
            node_sizes.append(max(400, int(200 * math.log2(hog_count + 1))))

        cycles = []
        try:
            cycles = nx.cycle_basis(bg)
        except Exception:
            pass
        cycle_nodes = set()
        for cyc in cycles:
            cycle_nodes.update(cyc)

        # Telomere blocks
        telomere_blocks = set()
        for cid, tels in self._child_telomeres.items():
            for tel_hog in tels:
                bid = self._hog_to_block.get(tel_hog)
                if bid:
                    telomere_blocks.add(bid)

        node_colors = []
        node_edge_colors = []
        for n in bg.nodes():
            if n in telomere_blocks:
                node_colors.append('#fff3cd')   # light yellow
                node_edge_colors.append('#e67e22')
            elif n in cycle_nodes:
                node_colors.append('#fce4ec')   # light pink
                node_edge_colors.append('#c62828')
            else:
                node_colors.append('#e3f2fd')   # light blue
                node_edge_colors.append('#1565c0')

        nx.draw_networkx_nodes(bg, pos, ax=ax, node_size=node_sizes,
                               node_color=node_colors, alpha=0.95,
                               edgecolors=node_edge_colors, linewidths=2.0)

        labels = {n: f'{n}\n({len(self._blocks.get(n, []))})' for n in bg.nodes()}
        nx.draw_networkx_labels(bg, pos, labels, ax=ax, font_size=10,
                                font_color='#1a1a2e', font_weight='bold')

        children = sorted(self.children())
        palette = ['#c0392b', '#2980b9', '#27ae60', '#f39c12',
                   '#8e44ad', '#16a085', '#d35400', '#2c3e50']
        child_colors = {c: palette[i % len(palette)]
                        for i, c in enumerate(children)}

        for h1, h2, data in bg.edges(data=True):
            colors = data.get('colors', set())
            child_ids = sorted(set(c for c, _ in colors))
            if len(child_ids) > 1:
                ax.plot([pos[h1][0], pos[h2][0]], [pos[h1][1], pos[h2][1]],
                        color='#1a1a2e', linewidth=3.5, alpha=0.9, zorder=1)
            elif len(child_ids) == 1:
                cid = child_ids[0]
                color = child_colors.get(cid, '#7f8c8d')
                ax.plot([pos[h1][0], pos[h2][0]], [pos[h1][1], pos[h2][1]],
                        color=color, linewidth=2.0, alpha=0.8, zorder=0)

        legend_handles = [
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='#fff3cd', markersize=12,
                       markeredgecolor='#e67e22', markeredgewidth=2,
                       label='telomere'),
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='#fce4ec', markersize=12,
                       markeredgecolor='#c62828', markeredgewidth=2,
                       label='in cycle'),
            plt.Line2D([0], [0], color='#1a1a2e', linewidth=3.5,
                       label='shared (multi-child)'),
        ]
        for cid in children:
            color = child_colors.get(cid, '#7f8c8d')
            legend_handles.append(
                plt.Line2D([0], [0], color=color, linewidth=2.5,
                           label=f'{cid} (unique)'))

        ax.legend(handles=legend_handles, loc='upper left', fontsize=8)
        ax.set_title(title or f'Block Graph: {self.hog_level}',
                     fontsize=12, fontweight='bold')
        ax.axis('off')
        fig.tight_layout()
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        logger.info("  [viz] block graph saved to %s", outpath)

    def draw_adjacency_heatmap(self, outpath: str, dpi: int = 200,
                               title: str = None):
        """Draw adjacency matrix heatmap: rows/cols = HOGs, cells = edge colors.

        Useful for seeing the global structure of shared vs unique adjacencies.
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.warning("matplotlib not available, skipping heatmap")
            return

        # Order HOGs by block, then by position within block
        if not hasattr(self, '_blocks'):
            self._build_synteny_blocks()

        hog_order = []
        for bid in sorted(self._blocks.keys()):
            hog_order.extend(self._blocks[bid])

        # Add HOGs not in any block
        all_hogs = set()
        for n in self.nodes():
            all_hogs.add(n)
        for h in sorted(all_hogs, key=str):
            if h not in set(hog_order):
                hog_order.append(h)

        n = len(hog_order)
        if n == 0 or n > 500:
            logger.info("  [viz] skipping heatmap: %d HOGs (too many or zero)", n)
            return

        hog_idx = {h: i for i, h in enumerate(hog_order)}

        # Build color matrix: 0=none, 1=unique, 2=shared
        mat = np.zeros((n, n), dtype=int)
        children = sorted(self.children())
        child_list = list(children)

        for h1, h2, data in self.edges(data=True):
            if h1 in hog_idx and h2 in hog_idx:
                i, j = hog_idx[h1], hog_idx[h2]
                colors = data.get('colors', set())
                if len(colors) > 1:
                    mat[i][j] = 2
                    mat[j][i] = 2
                else:
                    mat[i][j] = 1
                    mat[j][i] = 1

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cmap = matplotlib.colors.ListedColormap(['#ffffff', '#4ecdc4', '#ff6b6b'])
        ax.imshow(mat, cmap=cmap, interpolation='nearest', aspect='equal')

        # Block boundaries
        offset = 0
        for bid in sorted(self._blocks.keys()):
            blen = len(self._blocks[bid])
            ax.axhline(y=offset + blen - 0.5, color='gray', linewidth=0.5)
            ax.axvline(x=offset + blen - 0.5, color='gray', linewidth=0.5)
            offset += blen

        ax.set_title(title or f'Adjacency Matrix: {self.hog_level}',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('HOG index')
        ax.set_ylabel('HOG index')

        legend_handles = [
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#ffffff', markersize=10, label='no edge'),
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#4ecdc4', markersize=10, label='unique'),
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#ff6b6b', markersize=10, label='shared'),
        ]
        ax.legend(handles=legend_handles, loc='upper right', fontsize=8)

        fig.tight_layout()
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        logger.info("  [viz] adjacency heatmap saved to %s", outpath)

    def __repr__(self):
        return (f"ColoredGraph(hog_level={self.hog_level!r}, "
                f"nodes={self.node_count()}, edges={self.edge_count()}, "
                f"events={len(self.events)})")
