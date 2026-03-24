import os
import re
import pickle
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from constants import *
from matplotlib.patches import Patch

def sample_type(sample):
    name = SAMPLE_NAMES.get(sample, '')
    if '-' in name:
        suffix = name.split('-', 1)[1]
        return {'Silica': 'silica', 'Fresh': 'fresh', 'Herbarium': 'herbaria'}.get(suffix, 'unknown')
    return SAMPLE_TYPES.get(sample, 'unknown')


def display_name(sample):
    name = SAMPLE_NAMES.get(sample, sample)
    if '-' in name:
        genus, stype = name.split('-', 1)
        #bold and italicize the sample name
        return '$\\boldsymbol{' + genus + '}$-' + stype
    return name


CACHE = None


def parse_busco(path):
    content = open(path).read()
    lines   = content.split('\n')
    row     = next((l.strip() for l in lines if l.strip().startswith('C:') and '%' in l), None)
    if not row:
        return None
    matches = {k: re.search(rf'{k}:([\d.]+)%', row) for k in 'SDFM'}
    if not all(matches.values()):
        return None
    data = {k: float(matches[k].group(1)) for k in 'SDFM'}
    for line in lines:
        if 'Complete and single-copy BUSCOs (S)' in line:
            data['S_count'] = int(line.split()[0])
        elif 'Complete and duplicated BUSCOs (D)' in line:
            data['D_count'] = int(line.split()[0])
        elif 'Fragmented BUSCOs (F)' in line:
            data['F_count'] = int(line.split()[0])
        elif 'Missing BUSCOs (M)' in line:
            data['M_count'] = int(line.split()[0])
    return data


def parse_deamination(path):
    df = pd.read_csv(path, sep='\t')
    raw = {}
    for end in ('3p', '5p'):
        for std in ('+', '-'):
            sub = df[(df['End'] == end) & (df['Std'] == std)]
            for col in ('C>T', 'A>G', 'Total'):
                vals = sub[col].head(25).tolist() if col in sub.columns else []
                raw[(end, std, col)] = vals + [0] * (25 - len(vals))

    def combine(end, col):
        p, n   = raw[(end, '+', col)],     raw[(end, '-', col)]
        pt, nt = raw[(end, '+', 'Total')], raw[(end, '-', 'Total')]
        return [(a+b)/(c+d) if (c+d) > 0 else 0 for a, b, c, d in zip(p, n, pt, nt)]

    c3t = combine('3p', 'C>T'); c3t.reverse()
    a3g = combine('3p', 'A>G'); a3g.reverse()
    return {'5p_CtoT': combine('5p', 'C>T'), '5p_AtoG': combine('5p', 'A>G'),
            '3p_CtoT': c3t, '3p_AtoG': a3g}


def parse_inserts(path):
    df     = pd.read_csv(path)
    counts = df[df['insert_length'] >= 1]['insert_length'].value_counts().sort_index()
    x, y   = counts.index.values, counts.values
    # 0 is end to end and less than 0 is an overlap.
    # the data also fell off before 800 so no point in plotting all of it
    mask   = (x >= 1) & (x <= 800)
    return {'x': x[mask].tolist(), 'y': y[mask].tolist()}


def parse_transrate(path):
    df = pd.read_csv(path)
    return [float(df[c].iloc[0]) if c in df.columns else 0.0 for c in TRANSRATE_COLS]


def busco_path(sample):
    suffix = '_paired' if sample.startswith('WA') else ''
    name = f'{sample}{suffix}_busco'
    return os.path.join(DATA_DIR, 'busco', name, f'short_summary.specific.viridiplantae_odb10.{name}.txt')


def inserts_path(sample):
    suffix = '_paired' if sample.startswith('WA') else ''
    return os.path.join(DATA_DIR, 'inserts', f'{sample}{suffix}.csv')


# create a pickle of all the data needed. So that rerunning this for edits was faster than rescanning all the data
def build_cache():
    print("Building cache")
    cache = {}
    for sample in SAMPLES:
        entry = {}
        parsers = {
            'busco':       (busco_path(sample), parse_busco),
            'deamination': (os.path.join(DATA_DIR, 'deamination', sample, 'misincorporation.txt'), parse_deamination),
            'inserts':     (inserts_path(sample), parse_inserts),
            'transrate':   (os.path.join(DATA_DIR, 'transrate', f'{sample}.csv'), parse_transrate),
        }
        for key, (path, fn) in parsers.items():
            if os.path.exists(path):
                try:
                    entry[key] = fn(path)
                except Exception as e:
                    print(f"Warning: {key} parse failed for {sample}:\n{e}")
            else:
                print(f"Warning: {key} file not found for {sample}")
        cache[sample] = entry

    max_y = max((max(e['inserts']['y']) for e in cache.values()
                 if e.get('inserts') and e['inserts']['y']), default=0.0)
    cache['inserts_y'] = {'max_inserts_y': max_y * 1.1}

    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(cache, f)
    print(f"Cache written: {CACHE_FILE}")
    return cache


def load_cache():
    if os.path.exists(CACHE_FILE):
        print(f"Loading from cache")
        with open(CACHE_FILE, 'rb') as f:
            return pickle.load(f)
    return build_cache()


def get_cache():
    global CACHE
    if CACHE is None:
        CACHE = load_cache()
    return CACHE


def draw_busco(ax, sample, fonts, title=None, show_ylabel=True, **_):
    data   = get_cache().get(sample, {}).get('busco')
    pcts   = [data[k] for k in BUSCO_CATEGORIES]
    colors = [COLORS['busco'][k] for k in BUSCO_CATEGORIES]
    bars   = ax.bar(BUSCO_CATEGORIES, pcts, color=colors, alpha=0.8)

    for bar, pct in zip(bars, pcts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                f'{pct:.1f}%', ha='center', va='bottom',
                fontsize=fonts['bar_annotation'], fontweight='bold')

    ax.set_title(title or display_name(sample), fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks(range(len(BUSCO_CATEGORIES)))
    ax.set_xticklabels(BUSCO_LABELS, rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='x', labelsize=fonts['tick_label'], pad=8)
    ax.tick_params(axis='y', labelsize=fonts['tick_label'], pad=8)
    if show_ylabel:
        ax.set_ylabel('Percentage (%)', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y')


def draw_busco_cat(ax, category, category_label, fonts, show_ylabel=True, **_):
    cache    = get_cache()
    samples  = sorted(REPRESENTATIVE_SAMPLES, key=lambda s: TYPE_ORDER.index(sample_type(s)))
    vals     = [cache.get(s, {}).get('busco', {}).get(category, 0.0) or 0.0 for s in samples]
    colors   = [COLORS[sample_type(s)] for s in samples]
    xlabels  = [REP_SAMPLE_LABELS[s] for s in samples]
    bars     = ax.bar(range(len(samples)), vals, color=colors, alpha=0.8)

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                f'{val:.1f}%', ha='center', va='bottom',
                fontsize=fonts['bar_annotation'], fontweight='bold')

    ax.set_title(category_label, fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 110)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(xlabels, fontsize=fonts['tick_label'], rotation=45,
                       ha='right', rotation_mode='anchor')
    ax.tick_params(axis='y', labelsize=fonts['tick_label'])
    if show_ylabel:
        ax.set_ylabel('Percentage (%)', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y')


def draw_busco_all(ax, category, category_label, fonts, show_ylabel=True, **_):
    cache  = get_cache()
    groups = {t: [] for t in TYPE_ORDER}

    for s in SAMPLES:
        val = cache.get(s, {}).get('busco', {}).get(category)
        if val is not None:
            groups[sample_type(s)].append((s, val))
    for t in TYPE_ORDER:
        groups[t].sort(key=lambda x: x[1])

    x_pos, colors, values, labels = [], [], [], []
    x = 0
    for t in TYPE_ORDER:
        for s, val in groups[t]:
            x_pos.append(x)
            colors.append(COLORS[t])
            values.append(val)
            labels.append(SAMPLE_NAMES.get(s, s).split('-')[0])
            x += 1
        x += 1

    ax.bar(x_pos, values, color=colors, alpha=0.85, zorder=3)
    ax.set_title(category_label, fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 110)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=fonts['tick_label'], rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='y', labelsize=fonts['tick_label'])
    ax.tick_params(axis='x', pad=6)
    if show_ylabel:
        ax.set_ylabel('Percentage (%)', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y', zorder=0)


def draw_deamination(ax, sample, fonts, title=None, show_ylabel=True, show_legend=True, linewidth=2, **_):
    data = get_cache().get(sample, {}).get('deamination')
    x5   = list(range(1, 26))
    x3   = list(range(-25, 0))

    ax.plot(x5, data['5p_AtoG'], color='#56B4E9', label='A>I', linewidth=linewidth, alpha=0.8)
    ax.plot(x5, data['5p_CtoT'], color='#F04442', label='C>U', linewidth=linewidth, alpha=0.8)
    ax.plot(x3, data['3p_AtoG'], color='#56B4E9', linewidth=linewidth, alpha=0.8)
    ax.plot(x3, data['3p_CtoT'], color='#F04442', linewidth=linewidth, alpha=0.8)

    ax.set_title(title or display_name(sample), fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0.00, 0.002)
    ax.set_yticks([0.0000, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020])
    ax.set_xlim(-25.5, 25.5)
    ax.set_xticks([-25, -20, -15, -10, -5, -1, 1, 5, 10, 15, 20, 25])
    # 25 and -25 labels overlapped in the middle. so we had to get creative
    # this is also why you'll see sharey in the script
    ax.set_xticklabels(['1', '', '', '15', '', '', '', '', '-15', '', '', '-1'])
    ax.tick_params(axis='x', labelsize=fonts['tick_label'], pad=8)
    ax.tick_params(axis='y', labelsize=fonts['tick_label'], pad=8)
    ax.set_xlabel('Position', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    if show_ylabel:
        ax.set_ylabel('Misincorporation Frequency', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    if show_legend:
        ax.legend(fontsize=fonts['legend'], loc='upper right', framealpha=0.95, borderpad=1.0)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)


def draw_inserts(ax, sample, max_y, fonts, title=None, show_ylabel=True, linewidth=2, **_):
    data = get_cache().get(sample, {}).get('inserts')

    if data['x']:
        ax.plot(data['x'], data['y'], color='#000000', linewidth=linewidth)

    ax.set_title(title or display_name(sample), fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_xlim(1, 800)
    ax.set_ylim(-1, max_y)
    ax.set_xticks([1] + list(np.arange(100, 801, 100)))
    ax.tick_params(axis='x', rotation=45, labelsize=fonts['tick_label'], pad=8)
    ax.tick_params(axis='y', labelsize=fonts['tick_label'], pad=8)
    ax.set_xlabel('Insert Length (bp)', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    if show_ylabel:
        ax.set_ylabel('Count', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3)


def draw_transrate(ax, sample, fonts, title=None, show_ylabel=True, **_):
    scores = get_cache().get(sample, {}).get('transrate')
    x_pos  = np.arange(len(TRANSRATE_LABELS))
    bars   = ax.bar(x_pos, scores, color=[COLORS['transrate'][c] for c in TRANSRATE_COLS], alpha=0.8)

    for bar, score in zip(bars, scores):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{score:.3f}', ha='center', va='bottom',
                fontsize=fonts['bar_annotation'], fontweight='bold')

    ax.set_title(title or display_name(sample), fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 1.1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(TRANSRATE_LABELS, fontsize=fonts['tick_label'])
    ax.tick_params(axis='y', labelsize=fonts['tick_label'])
    if show_ylabel:
        ax.set_ylabel('Score', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y')


def draw_transrate_score(ax, score_col, score_label, fonts, show_ylabel=True, **_):
    score_idx = TRANSRATE_COLS.index(score_col)
    cache     = get_cache()
    samples   = sorted(REPRESENTATIVE_SAMPLES, key=lambda s: TYPE_ORDER.index(sample_type(s)))
    vals      = [cache.get(s, {}).get('transrate', [0]*4)[score_idx] for s in samples]
    colors    = [COLORS[sample_type(s)] for s in samples]
    xlabels   = [REP_SAMPLE_LABELS[s] for s in samples]
    bars      = ax.bar(range(len(samples)), vals, color=colors, alpha=0.8)

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', va='bottom',
                fontsize=fonts['bar_annotation'], fontweight='bold')

    ax.set_title(score_label, fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 1.1)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(xlabels, fontsize=fonts['tick_label'], rotation=45,
                       ha='right', rotation_mode='anchor')
    ax.tick_params(axis='y', labelsize=fonts['tick_label'])
    if show_ylabel:
        ax.set_ylabel('Score', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y')


def draw_transrate_all(ax, score_col, score_label, fonts, show_ylabel=True):
    score_idx = TRANSRATE_COLS.index(score_col)
    cache     = get_cache()
    groups    = {t: [] for t in TYPE_ORDER}

    for s in SAMPLES:
        scores = cache.get(s, {}).get('transrate')
        if scores:
            groups[sample_type(s)].append((s, scores[score_idx]))
    for t in TYPE_ORDER:
        groups[t].sort(key=lambda x: x[1])

    x_pos  = []
    colors = []
    values = []
    labels = []
    x      = 0

    for t in TYPE_ORDER:
        for s, val in groups[t]:
            x_pos.append(x)
            colors.append(COLORS[t])
            values.append(val)
            labels.append(SAMPLE_NAMES.get(s, s).split('-')[0])
            x += 1
        x += 1

    ax.bar(x_pos, values, color=colors, alpha=0.85, zorder=3)
    ax.set_title(score_label, fontsize=fonts['title'], fontweight='bold', pad=fonts['title_pad'])
    ax.set_ylim(0, 1.1)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=fonts['tick_label'], rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='y', labelsize=fonts['tick_label'])
    ax.tick_params(axis='x', pad=6)
    if show_ylabel:
        ax.set_ylabel('Score', fontsize=fonts['axis_label'], labelpad=fonts['labelpad'])
    ax.grid(True, alpha=0.3, axis='y', zorder=0)


# in the concat rep image, busco and transrate are 4 columns instead of 3
def draw_4col(fig, gs, row_idx, items, draw_fn, fonts):
    col_slices = [slice(0, 3), slice(3, 6), slice(6, 9), slice(9, 12)]
    first_ax   = None
    for col_idx, (key, label) in enumerate(items):
        ax = fig.add_subplot(gs[row_idx, col_slices[col_idx]])
        draw_fn(ax, key, label, fonts, show_ylabel=(col_idx == 0))
        if col_idx > 0:
            ax.set_ylabel('')
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            first_ax = ax
    return first_ax


# the 4 plots of the 3 rep samples
def plot_rep(plot_type, output_dir):
    print(f"{plot_type.capitalize():<12} - representative")
    sub_dir = os.path.join(output_dir, plot_type)
    os.makedirs(sub_dir, exist_ok=True)

    max_y = get_cache()['inserts_y']['max_inserts_y'] if plot_type == 'inserts' else None
    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE_REP, sharey=True)

    lw = 8 if plot_type in ('deamination', 'inserts') else 2
    draw_fn = {
        'busco'      : lambda ax, s, i: draw_busco(ax, s, FONTS_CONCAT, show_ylabel=(i==0)),
        'deamination': lambda ax, s, i: draw_deamination(ax, s, FONTS_CONCAT, show_ylabel=(i==0), show_legend=(i==0), linewidth=lw),
        'inserts'    : lambda ax, s, i: draw_inserts(ax, s, max_y, FONTS_CONCAT, show_ylabel=(i==0), linewidth=lw),
        'transrate'  : lambda ax, s, i: draw_transrate(ax, s, FONTS_CONCAT, show_ylabel=(i==0)),
    }[plot_type]

    for i, sample in enumerate(REPRESENTATIVE_SAMPLES):
        draw_fn(axes[i], sample, i)
        if i > 0:
            axes[i].set_ylabel('')
            plt.setp(axes[1].get_yticklabels(), visible=False)

    plt.tight_layout(pad=8.0, h_pad=4.0, w_pad=3.0)
    plt.savefig(os.path.join(sub_dir, f'{plot_type}_representative.png'), dpi=300, bbox_inches='tight')
    plt.close()

# transrate needed its own function because it has one that is all the individuals combined
# and this is a seperate one thats grouped by score type. snuc scov etc 
def plot_transrate_scores(output_dir):

    print(f"{'Transrate':<12} - scores")
    sub_dir = os.path.join(output_dir, 'transrate')
    os.makedirs(sub_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_TRANSRATE_SCORES)
    for ax, (col, label) in zip(axes.flatten(), zip(TRANSRATE_COLS, TRANSRATE_LABELS)):
        draw_transrate_all(ax, col, label, FONTS_IND)

    legend_patches = [
        Patch(facecolor=COLORS['fresh'],    alpha=0.85, label='Fresh'),
        Patch(facecolor=COLORS['silica'],   alpha=0.85, label='Silica'),
        Patch(facecolor=COLORS['herbaria'], alpha=0.85, label='Herbarium'),
    ]
    fig.legend(handles=legend_patches, fontsize=FONTS_IND['legend'],
               loc='upper center', ncol=3, framealpha=0.8, bbox_to_anchor=(0.5, 1.01))

    plt.tight_layout(pad=2.5)
    plt.savefig(os.path.join(sub_dir, 'transrate_scores.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_busco_categories(output_dir):

    print(f"{'Busco':<12} - categories")
    sub_dir = os.path.join(output_dir, 'busco')
    os.makedirs(sub_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_BUSCO_SCORES)
    for ax, (cat, label) in zip(axes.flatten(), zip(BUSCO_CATEGORIES, BUSCO_LABELS)):
        draw_busco_all(ax, cat, label.replace('\n', ' '), FONTS_IND)

    legend_patches = [
        Patch(facecolor=COLORS['fresh'],    alpha=0.85, label='Fresh'),
        Patch(facecolor=COLORS['silica'],   alpha=0.85, label='Silica'),
        Patch(facecolor=COLORS['herbaria'], alpha=0.85, label='Herbarium'),
    ]
    fig.legend(handles=legend_patches, fontsize=FONTS_IND['legend'],
               loc='upper center', ncol=3, framealpha=0.8, bbox_to_anchor=(0.5, 1.01))

    plt.tight_layout(pad=2.5)
    plt.savefig(os.path.join(sub_dir, 'busco_categories.png'), dpi=300, bbox_inches='tight')
    plt.close()


# this is figure 2 (i think 2?) with all of the rep plots in 4 rows.
def concat_rep(output_dir, plot_types):
    print("Concatenated Representative Plots")
    ordered = [pt for pt in CONCAT_PLOT_ORDER if pt in plot_types]
    if len(ordered) < 2:
        return

    max_y          = get_cache()['inserts_y']['max_inserts_y'] if 'inserts' in ordered else None
    fig            = plt.figure(figsize=FIGSIZE_CONCAT)
    gs             = gridspec.GridSpec(len(ordered), 12, figure=fig)
    col3           = [slice(0, 4), slice(4, 8), slice(8, 12)]
    row_first_axes = []

    for row_idx, plot_type in enumerate(ordered):
        if plot_type == 'transrate':
            first = draw_4col(fig, gs, row_idx, list(zip(TRANSRATE_COLS, TRANSRATE_LABELS)), draw_transrate_score, FONTS_CONCAT)
            row_first_axes.append(first)

        elif plot_type == 'busco':
            first = draw_4col(fig, gs, row_idx, list(zip(BUSCO_CATEGORIES, BUSCO_LABELS)), draw_busco_cat, FONTS_CONCAT)
            row_first_axes.append(first)

        else:
            lw = {'linewidth': 8} if plot_type in ('deamination', 'inserts') else {}
            for col_idx, sample in enumerate(CONCAT_SAMPLES):
                ax = fig.add_subplot(gs[row_idx, col3[col_idx]])
                if plot_type == 'deamination':
                    draw_deamination(ax, sample, FONTS_CONCAT, show_ylabel=(col_idx == 0), show_legend=(col_idx == 0), **lw)
                elif plot_type == 'inserts':
                    draw_inserts(ax, sample, max_y, FONTS_CONCAT, show_ylabel=(col_idx == 0), **lw)
                if col_idx > 0:
                    ax.set_ylabel('')
                    plt.setp(ax.get_yticklabels(), visible=False)
                else:
                    row_first_axes.append(ax)

    plt.tight_layout(pad=8.0, h_pad=4.0, w_pad=3.0, rect=[0.04, 0.08, 0.92, 0.96])

    for row_idx, plot_type in enumerate(ordered):
        pos = row_first_axes[row_idx].get_position()
        fig.text(0.03, pos.y1 + 0.005, chr(65 + row_idx), fontsize=FONT_ROW_LABEL, fontweight='bold', ha='center', va='top')
        fig.text(pos.x0 - 0.09, pos.y0 + (pos.y1 - pos.y0) / 2, ROW_TITLES[plot_type], fontsize=FONT_ROW_TITLE, fontweight='bold', ha='center', va='center', rotation=90)

    plt.savefig(os.path.join(output_dir, 'representative_plots_concatenated.svg'), format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'representative_plots_concatenated.png'), format='png', bbox_inches='tight')
    plt.close()


# individual thread/worker for individual plots, based on how many threads wanted
def individual_worker(args):
    plot_type, sample, output_dir, max_y = args
    cache   = get_cache()
    fig, ax = plt.subplots(1, 1, figsize=FIGSIZE[plot_type])
    title   = f'{display_name(sample)}{TITLE_SUFFIX[plot_type]}'

    if plot_type == 'busco':
        draw_busco(ax, sample, FONTS_IND, title=title)
    elif plot_type == 'deamination':
        draw_deamination(ax, sample, FONTS_IND, title=title)
    elif plot_type == 'inserts':
        draw_inserts(ax, sample, max_y, FONTS_IND, title=title)
    elif plot_type == 'transrate':
        draw_transrate(ax, sample, FONTS_IND, title=title)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, plot_type, 'individual', f'{plot_type}_{sample}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    return


# makes a powerpoint slide of all the individual plots concatenated together
def plot_all_samples(plot_type, output_dir):
    print(f"{plot_type.capitalize():<12} - all samples")
    max_y   = get_cache()['inserts_y']['max_inserts_y'] if plot_type == 'inserts' else None
    samples = sorted(SAMPLES.keys(),
                     key=lambda s: (TYPE_ORDER.index(sample_type(s)), s))

    fig, axes = plt.subplots(SLIDE_NROWS, SLIDE_NCOLS, figsize=SLIDE_FIGSIZE)
    axes_flat = axes.flatten()

    for i, (ax, sample) in enumerate(zip(axes_flat, samples)):
        show_ylabel = (i % SLIDE_NCOLS == 0)
        if plot_type == 'busco':
            draw_busco(ax, sample, FONTS_SLIDE, show_ylabel=show_ylabel)
        elif plot_type == 'deamination':
            draw_deamination(ax, sample, FONTS_SLIDE,
                             show_ylabel=show_ylabel, show_legend=False)
        elif plot_type == 'inserts':
            draw_inserts(ax, sample, max_y, FONTS_SLIDE, show_ylabel=show_ylabel)
        elif plot_type == 'transrate':
            draw_transrate(ax, sample, FONTS_SLIDE, show_ylabel=show_ylabel)

    for ax in axes_flat[len(samples):]:
        ax.set_visible(False)

    plt.tight_layout(pad=0.8, h_pad=0.6, w_pad=0.4)
    out_path = os.path.join(output_dir, plot_type, f'{plot_type}_all_samples.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

# individual plots for each sample of each analysis, transrate, busco etc
def plot_individual(plot_type, output_dir, threads=16):
    print(f"{plot_type.capitalize():<12} - individual")
    os.makedirs(os.path.join(output_dir, plot_type, 'individual'), exist_ok=True)
    max_y = get_cache()['inserts_y']['max_inserts_y'] if plot_type == 'inserts' else None
    args  = [(plot_type, s, output_dir, max_y) for s in SAMPLES]
    with mp.Pool(processes=threads) as pool:
        pool.map(individual_worker, args)
    plot_all_samples(plot_type, output_dir)
