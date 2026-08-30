#!/bin/env python
# coding: utf-8
import re
"""Global style parameters for vis plotting. Analogous to matplotlib rcParams."""

STYLE = {
    'fontsize': {
        'chrom_label': 1.0,      # multiplier on base fontsize
        'axis_label': 1.5,
        'tick_label': 0.6,
        'title': 1.5,
        'cbar_label': 0.8,
        'cbar_tick': 0.5,
    },
    'linewidth': {
        'chr_boundary': 1.0,
        'chr_inner': 0.8,
        'sym_triangle': 1.2,
    },
    'color': {
        'chr_boundary': '#444444',       # deep gray
        'chr_inner': '#CCCCCC',          # light gray
        'symbols': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3'],
    },
    'dashes': {
        'chr_inner': (0, (15, 12)),       # custom dash for inner lines
    },
    'symbol': {
        'alpha': 0.7,
        'height_ratio': 0.02,   # fraction of total axis length
        'width_ratio': 0.01,
        'jitter_threshold': 0.5,  # fraction of width to trigger jitter
        'jitter_step': 0.8,       # fraction of width per step
        'edge_gap_ratio': 0.15,   # gap between triangle tip and plot edge (fraction of h)
    },
}


def format_species_label(name):
	"""Latin-name pattern (Genus_species, all-lower tail) -> (text with
	underscores as spaces, 'italic'); else (name, 'normal')."""
	if re.match(r'^[A-Z][a-z]+(_[a-z]+)+$', name):
		return (name.replace('_', ' '), 'italic')
	return (name, 'normal')


def style_species(t, keep_style=False, **kw):
	"""Style a matplotlib Text in place by the species-name pattern.
	t: Text instance (from set_xlabel/set_ylabel/get_yticklabels/
	legend.get_texts(), ...). Latin-name text -> italic + underscores
	as spaces; otherwise left as-is when keep_style=True, else
	fontstyle='normal'. Extra attributes go to t.set(**kw).
	Returns t (chainable)."""
	text, style = format_species_label(t.get_text())
	if text != t.get_text():
		t.set_text(text)
	if not keep_style:
		t.set_fontstyle(style)
	if kw:
		t.set(**kw)
	return t


def set_species_ticklabels(ax, names, axis='y', **kw):
	"""One-call species tick labels: set ticks + labels + per-label
	italic/normal via style_species. MUST be called after operations
	that rebuild tick labels (yaxis.tick_right(), tick_params with
	labelleft changes, ...). kw forwards to set_*ticklabels
	(e.g. fontsize)."""
	getter = ax.get_yticklabels if axis == 'y' else ax.get_xticklabels
	setter = ax.set_yticklabels if axis == 'y' else ax.set_xticklabels
	ticksetter = ax.set_yticks if axis == 'y' else ax.set_xticks
	ticksetter(range(len(names)))
	labels = [format_species_label(n) for n in names]
	setter([t for t, s in labels], **kw)
	for lb, (t, s) in zip(getter(), labels):
		lb.set_fontstyle(s)
