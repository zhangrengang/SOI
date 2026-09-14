#!/bin/env python
# coding: utf-8
import re
"""Global style parameters for vis plotting. Analogous to matplotlib rcParams.
Importing this module applies the global matplotlib settings (Agg backend,
TrueType output, Arial-first sans-serif)."""


DEFAULT_FONT = 'Arial'


def setup_mpl():
	import matplotlib as mpl
	mpl.use('Agg')
	mpl.rcParams.update({
		'pdf.fonttype': 42,
		'ps.fonttype': 42,
		'font.family': 'sans-serif',
		'font.sans-serif': [DEFAULT_FONT, 'DejaVu Sans'],  # DEFAULT_FONT optional; DejaVu fallback
	})


setup_mpl()

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


# Latin binomial: Genus_species
PATTERN_SPECIES = r'^[A-Z][a-z]+(_[a-z]+)+$'
# Latin binomial with an optional trailing sample suffix (e.g.
# Helicteres_baruensis_1, Reevesia_pubescens_WGS, ..._RNA)
PATTERN_SPECIES_SUFFIX = r'^[A-Z][a-z]+(_[a-z]+)+(_[A-Za-z0-9]+)?$'


def format_species_label(name, suffix=True):
	"""Latin-name pattern (Genus_species, all-lower tail) -> (text with
	underscores as spaces, 'italic'); else (name, 'normal').

	suffix=True (default): also italicise binomials carrying a trailing
	sample suffix (_1 / _WGS / _RNA ...), which is how multiple datasets
	of one species are labelled.
	"""
	if re.match(PATTERN_SPECIES_SUFFIX if suffix else PATTERN_SPECIES, name):
		return (name.replace('_', ' '), 'italic')
	return (name, 'normal')


def style_species(t, keep_style=False, suffix=True, **kw):
	"""Style a matplotlib Text in place by the species-name pattern.
	t: Text instance (from set_xlabel/set_ylabel/get_yticklabels/
	legend.get_texts(), ...). Latin-name text -> italic + underscores
	as spaces; otherwise left as-is when keep_style=True, else
	fontstyle='normal'. Extra attributes go to t.set(**kw).
	Returns t (chainable)."""
	text, style = format_species_label(t.get_text(), suffix=suffix)
	if text != t.get_text():
		t.set_text(text)
	if not keep_style:
		t.set_fontstyle(style)
	if kw:
		t.set(**kw)
	return t


def set_species_ticklabels(ax, names, axis='y', suffix=True, **kw):
	"""One-call species tick labels: set ticks + labels + per-label
	italic/normal via style_species. MUST be called after operations
	that rebuild tick labels (yaxis.tick_right(), tick_params with
	labelleft changes, ...). kw forwards to set_*ticklabels
	(e.g. fontsize)."""
	getter = ax.get_yticklabels if axis == 'y' else ax.get_xticklabels
	setter = ax.set_yticklabels if axis == 'y' else ax.set_xticklabels
	ticksetter = ax.set_yticks if axis == 'y' else ax.set_xticks
	ticksetter(range(len(names)))
	labels = [format_species_label(n, suffix=suffix) for n in names]
	setter([t for t, s in labels], **kw)
	for lb, (t, s) in zip(getter(), labels):
		lb.set_fontstyle(s)


def format_species_title(s, sep=' vs ', suffix=True):
	"""'Genus_a vs Genus_b' style titles: if BOTH sides match the Latin
	pattern, return (text with underscores as spaces, 'italic');
	otherwise fall back to single-name format_species_label (handles
	bare species names and node names like N4)."""
	parts = s.split(sep)
	if len(parts) == 2:
		fmts = [format_species_label(p, suffix=suffix) for p in parts]
		if all(st == 'italic' for t, st in fmts):
			return (sep.join(t for t, st in fmts), 'italic')
	return format_species_label(s, suffix=suffix)
