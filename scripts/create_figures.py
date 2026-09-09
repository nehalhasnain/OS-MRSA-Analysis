#!/usr/bin/env python3
"""
OS-MRSA Project: Publication Figures Reproduction Script
Target Journal: One Health Advances (BioMed Central / Springer Nature)

Description:
  Reproduces all main and supplementary figures (Figures 1–4, Figures S1–S4)
  from the de-identified minimal dataset.

Requirements:
  Python >= 3.8
  numpy, pandas, matplotlib, scipy

Usage:
  python create_figures.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle, FancyArrow, Patch
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import beta

# Determine project directories
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

# Output directory
OUT = PROJECT_ROOT / 'figures'
MAIN = OUT / 'Main_Figures'
SUPP = OUT / 'Supplementary_Figures'
DATA_OUT = OUT / 'Source_Data'
for d in (MAIN, SUPP, DATA_OUT):
    d.mkdir(parents=True, exist_ok=True)

# Locate input dataset
possible_paths = [
    PROJECT_ROOT / 'data' / 'Additional_file_1_Deidentified_Minimal_Dataset.csv',
    SCRIPT_DIR / 'Additional_file_1_Deidentified_Minimal_Dataset.csv',
    Path('data/Additional_file_1_Deidentified_Minimal_Dataset.csv'),
    Path('Additional_file_1_Deidentified_Minimal_Dataset.csv')
]

data_file = None
for p in possible_paths:
    if p.exists():
        data_file = p
        break

if not data_file:
    raise FileNotFoundError(
        "Could not find 'Additional_file_1_Deidentified_Minimal_Dataset.csv'. "
        "Ensure it exists in the 'data/' directory."
    )

print(f"Loading dataset: {data_file}")

# Typography configuration (fallback gracefully if system lacks local TTF paths)
font_candidates = [
    '/System/Library/Fonts/Supplemental/Times New Roman.ttf',
    '/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf',
    '/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf',
    '/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf',
    'C:\\Windows\\Fonts\\times.ttf'
]
for fc in font_candidates:
    if Path(fc).exists():
        try:
            fm.fontManager.addfont(fc)
        except Exception:
            pass

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Liberation Serif'],
    'font.size': 8.5,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'axes.linewidth': 0.65,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

BLUE, AMBER, TEAL, GREY, RED = '#0072B2', '#E69F00', '#009E73', '#747474', '#C63D2F'
LIGHT, MID = '#E5E5E5', '#4D4D4D'
SUS, INT, RES = BLUE, '#F2A93B', '#D84315'

def save(fig, path):
    fig.savefig(path.with_suffix('.png'), dpi=300, bbox_inches='tight', pad_inches=0.04, facecolor='white')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight', pad_inches=0.04, facecolor='white')
    try:
        fig.savefig(path.with_suffix('.tiff'), dpi=600, bbox_inches='tight', pad_inches=0.04, facecolor='white')
    except Exception:
        pass
    plt.close(fig)

def cp(k, n):
    lo = 0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1 if k == n else beta.ppf(0.975, k + 1, n - k)
    return np.array([lo, hi]) * 100

def stat_text(k, n):
    return f'{k}/{n}'

def clean(x):
    return str(x).strip() if pd.notna(x) else np.nan

df = pd.read_csv(data_file, encoding='utf-8-sig')
for c in ['S_aureus_Confirmed', 'PCR_mecA', 'Cefoxitin_FOX', 'Health_status', 'Prior_Antibiotics_Use']:
    if c in df.columns:
        df[c] = df[c].map(clean)

sa = df[df.S_aureus_Confirmed.eq('Yes')].copy()
assert len(sa) == 24, f"Expected 24 S. aureus isolates, found {len(sa)}"

# Derived resistance burden
sa['Resistance_classes'] = sa['Resistance_class_count'].astype(int)
burden = sa[['Sample_ID', 'PCR_mecA', 'Resistance_classes']].copy()
burden['mecA_status'] = burden['PCR_mecA'].map({'Yes': 'mecA-positive', 'No': 'mecA-negative'})

# ==============================================================================
# FIGURE 1: Compact prevalence forest plot by animal species
# ==============================================================================
print("Generating Figure 1: Species prevalence forest plot...")
species = ['Cat', 'Dog', 'Horse']
outcomes = [(r'$\it{S.\ aureus}$', BLUE, 'S_aureus_Confirmed'), (r'$\it{mecA}$-positive MRSA', RED, 'PCR_mecA')]
rows = []
for sp in species:
    d = df[df.Species.eq(sp)]
    for label, color, col in outcomes:
        k = int(d[col].eq('Yes').sum())
        n = len(d)
        lo, hi = cp(k, n)
        rows.append((sp, label, k, n, 100 * k / n, lo, hi))

f1 = pd.DataFrame(rows, columns=['species', 'outcome', 'positive', 'denominator', 'percent', 'lower95', 'upper95'])
f1.to_csv(DATA_OUT / 'Figure_1_species_prevalence.csv', index=False)

fig, ax = plt.subplots(figsize=(6.5, 3.55))
ybase = np.array([2, 1, 0])
offsets = [0.14, -0.14]
for off, (label, color, _) in zip(offsets, outcomes):
    s = f1[f1.outcome.eq(label)].set_index('species').loc[species]
    yy = ybase + off
    xx = s.percent.to_numpy()
    ci = np.vstack([xx - s.lower95.to_numpy(), s.upper95.to_numpy() - xx])
    ax.errorbar(xx, yy, xerr=ci, fmt='o', color=color, ecolor=color, markersize=5.5, capsize=2.4, lw=1.05, label=label, zorder=3)
    for x, y, k, n in zip(xx, yy, s.positive, s.denominator):
        ax.text(106, y, stat_text(k, n), va='center', fontsize=8)

ax.axvline(0, color=MID, lw=0.6)
ax.set_xlim(0, 123)
ax.set_ylim(-0.65, 2.65)
ax.set_yticks(ybase, species)
ax.set_xticks(range(0, 101, 20), [f'{x}%' for x in range(0, 101, 20)])
ax.set_xlabel('Prevalence (exact 95% confidence interval)')
ax.grid(axis='x', color=LIGHT, lw=0.55)
ax.set_axisbelow(True)
ax.text(0, 1.04, 'Prevalence by animal species', transform=ax.transAxes, fontweight='bold', fontsize=10.5)
ax.text(106, 2.47, 'Events / total', fontsize=8, fontweight='bold')
ax.text(0, -0.22, 'Horse denominator was 1; estimate is descriptive.', transform=ax.transAxes, fontsize=7.2, color=MID)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.legend(frameon=False, loc='lower left', ncol=2, bbox_to_anchor=(0, 1.05), handletextpad=0.4, columnspacing=1.3)
fig.tight_layout()
save(fig, MAIN / 'Figure_1_species_prevalence')

# ==============================================================================
# FIGURE 2: Subgroup prevalence forest plot (Health status & Prior antibiotics)
# ==============================================================================
print("Generating Figure 2: Subgroup prevalence forest plot...")
specs = [
    ('Health status', 'Health_status', [('Apparently healthy', 'Apparently Healthy'), ('Sick', 'Sick')]),
    ('Prior antibiotic use', 'Prior_Antibiotics_Use', [('Yes', 'Yes'), ('No', 'No'), ('Unknown', 'Unknown')])
]
allrows = []
for header, col, groups in specs:
    for dis, value in groups:
        d = df[df[col].eq(value)]
        for label, color, pcr in outcomes:
            k = int(d[pcr].eq('Yes').sum())
            n = len(d)
            lo, hi = cp(k, n)
            allrows.append((header, dis, label, k, n, 100 * k / n, lo, hi))

f2 = pd.DataFrame(allrows, columns=['panel', 'group', 'outcome', 'positive', 'denominator', 'percent', 'lower95', 'upper95'])
f2.to_csv(DATA_OUT / 'Figure_2_subgroup_prevalence.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.7), sharex=True)
for ax, (header, col, groups) in zip(axes, specs):
    order = [g[0] for g in groups]
    yy = np.arange(len(order))[::-1]
    for off, (label, color, _) in zip([0.12, -0.12], outcomes):
        s = f2[(f2.panel.eq(header)) & (f2.outcome.eq(label))].set_index('group').loc[order]
        xx = s.percent.to_numpy()
        ci = np.vstack([xx - s.lower95.to_numpy(), s.upper95.to_numpy() - xx])
        ax.errorbar(xx, yy + off, xerr=ci, fmt='o', color=color, ecolor=color, markersize=5, capsize=2.2, lw=1.0, zorder=3)
        for x, y, k, n in zip(xx, yy + off, s.positive, s.denominator):
            ax.text(min(x + 2.0, 98), y, stat_text(k, n), va='center', fontsize=7.5)
    ax.set_title(header, loc='left', fontweight='bold', fontsize=10, pad=7)
    ax.set(yticks=yy, yticklabels=order, xlim=(0, 105), ylim=(-0.6, len(order) - 0.35))
    ax.set_xticks(range(0, 101, 25), [f'{x}%' for x in range(0, 101, 25)])
    ax.grid(axis='x', color=LIGHT, lw=0.55)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0)

axes[0].set_xlabel('Prevalence (95% CI)')
axes[1].set_xlabel('Prevalence (95% CI)')
axes[0].legend(
    handles=[
        plt.Line2D([], [], marker='o', color=BLUE, linestyle='None', label=r'$\it{S.\ aureus}$'),
        plt.Line2D([], [], marker='o', color=RED, linestyle='None', label=r'$\it{mecA}$-positive MRSA')
    ],
    frameon=False, ncol=2, loc='lower left', bbox_to_anchor=(0, 1.17), columnspacing=1.2
)
fig.text(0.01, 0.01, 'Exact binomial confidence intervals. Unknown prior antibiotic use: n = 5.', fontsize=7.1, color=MID)
fig.subplots_adjust(top=0.76, bottom=0.20, wspace=0.36, left=0.10, right=0.98)
save(fig, MAIN / 'Figure_2_subgroup_prevalence_CI')

# ==============================================================================
# FIGURE 3: Genotype-phenotype discordance matrix
# ==============================================================================
print("Generating Figure 3: Genotype-phenotype discordance matrix...")
matrix = pd.DataFrame({
    'mecA_status': ['positive', 'positive', 'negative', 'negative'],
    'cefoxitin_status': ['resistant', 'susceptible', 'resistant', 'susceptible'],
    'classification': ['Typical MRSA', 'OS-MRSA', 'Cefoxitin-resistant mecA-negative', 'MSSA'],
    'count': [12, 6, 3, 3],
    'percent_of_24': [50.0, 25.0, 12.5, 12.5]
})
matrix.to_csv(DATA_OUT / 'Figure_3_genotype_phenotype.csv', index=False)

fig, ax = plt.subplots(figsize=(7.0, 4.5))
values = [[12, 6], [3, 3]]
names = [['Typical MRSA', 'OS-MRSA'], ['Cefoxitin-resistant\nmecA-negative', 'MSSA']]
colors = [[BLUE, AMBER], [GREY, TEAL]]
for i in range(2):
    for j in range(2):
        ax.add_patch(Rectangle((j, 1 - i), 1, 1, fc=colors[i][j], ec='white', lw=3))
        tc = 'black' if (i, j) == (0, 1) else 'white'
        v = values[i][j]
        ax.text(j + 0.5, 1 - i + 0.62, f'{v}/24', ha='center', va='center', fontsize=18, fontweight='bold', color=tc)
        ax.text(j + 0.5, 1 - i + 0.42, f'({100 * v / 24:.1f}%)', ha='center', va='center', fontsize=9.5, color=tc)
        ax.text(j + 0.5, 1 - i + 0.18, names[i][j], ha='center', va='center', fontsize=9, color=tc)

ax.set(xlim=(0, 2), ylim=(-0.29, 2), aspect='equal')
ax.set_xticks([0.5, 1.5], ['Cefoxitin resistant', 'Cefoxitin susceptible'])
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xlabel('Disc-diffusion interpretation', labelpad=8)
ax.set_yticks([1.5, 0.5], [r'$\it{mecA}$ positive (n = 18)', r'$\it{mecA}$ negative (n = 6)'])
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

ax.text(1, -0.18, r'OS-MRSA: 6 of 18 (33.3%) $\it{mecA}$-positive isolates were cefoxitin susceptible.', ha='center', va='center', fontsize=9.5, fontweight='bold')
fig.subplots_adjust(left=0.24, right=0.97, top=0.83, bottom=0.15)
save(fig, MAIN / 'Figure_3_genotype_phenotype_discordance')

# ==============================================================================
# FIGURE 4: Ordered locus map on NODE_13
# ==============================================================================
print("Generating Figure 4: Ordered locus map on NODE_13...")
features = [
    ('ccrA', 3736, 5085, '+', TEAL, 'ccr'),
    ('ccrB', 5107, 6735, '+', TEAL, 'ccr'),
    ('ISSep1', 8662, 10185, '+', '#8F56A6', 'is'),
    ('ΔmecR1', 10408, 11394, '-', AMBER, 'mec'),
    ('mecA', 11494, 13500, '+', RED, 'mec'),
    ('ugpQ', 14071, 14814, '-', '#4D95C6', 'other'),
    ('IS431mec', 16036, 16710, '+', '#8F56A6', 'is'),
    ('hypothetical ORF', 17396, 18652, '+', '#C4C4C4', 'other'),
    ('rlmH', 18974, 19453, '-', '#F28E1C', 'other')
]
pd.DataFrame(features, columns=['feature', 'start_bp', 'end_bp', 'strand', 'color', 'category']).to_csv(DATA_OUT / 'Figure_4_NODE13_locus_coordinates.csv', index=False)

fig, ax = plt.subplots(figsize=(8.1, 2.9))
base, end = 3300, 20200
ax.plot([0, (end - base) / 1000], [0.5, 0.5], color=MID, lw=0.65, zorder=0)
for i, (name, start, stop, strand, color, cat) in enumerate(features):
    x = (start - base) / 1000
    w = (stop - start) / 1000
    dx = w if strand == '+' else -w
    x0 = x if strand == '+' else x + w
    ax.add_patch(FancyArrow(x0, 0.5, dx, 0, width=0.17, head_width=0.17, head_length=min(0.20, w * 0.28), length_includes_head=True, fc=color, ec='#222222', lw=0.4))
    ly = 0.83 if name in ['ccrA', 'ISSep1', 'mecA', 'IS431mec', 'rlmH'] else 0.18
    ax.text(x + w / 2, ly, name, ha='center', va='center', fontsize=7.6, fontstyle='italic' if name in ['ccrA', 'ccrB', 'ΔmecR1', 'mecA', 'ugpQ', 'rlmH'] else 'normal')

# Context brackets
x1 = (10408 - base) / 1000
x2 = (13500 - base) / 1000
ax.plot([x1, x1, x2, x2], [0.98, 1.04, 1.04, 0.98], color=MID, lw=0.65)
ax.text((x1 + x2) / 2, 1.13, 'Class B mec complex', ha='center', fontsize=8)

# Scale bar
ax.plot([1, 3], [-0.12, -0.12], color='black', lw=1.4)
ax.text(2, -0.24, '2 kb', ha='center', fontsize=7.5)
ax.set(xlim=(0, (end - base) / 1000), ylim=(-0.33, 1.28), yticks=[], xlabel='Position on NODE_13 (kb)')
ticks = np.arange(4, 21, 2)
ax.set_xticks(ticks - base / 1000, ticks)
ax.tick_params(axis='x', length=3)
ax.spines[['left', 'right', 'top']].set_visible(False)
ax.text(0, 1.32, 'SCCmec-associated locus on NODE_13', fontsize=10.5, fontweight='bold', transform=ax.transAxes)
handles = [Patch(fc=TEAL, label='ccr genes'), Patch(fc=RED, label='mec complex'), Patch(fc='#8F56A6', label='insertion sequence'), Patch(fc='#C4C4C4', label='other ORF')]
ax.legend(handles=handles, frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.58), handlelength=1, columnspacing=1.1)
fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.28)
save(fig, MAIN / 'Figure_4_SCCmec_associated_locus_NODE13')

# ==============================================================================
# FIGURE S1: AST profile stacked horizontal bars
# ==============================================================================
print("Generating Figure S1: AST profile...")
agents = [
    ('Erythromycin_E', 'Erythromycin'), ('Ciprofloxacin_CIP', 'Ciprofloxacin'),
    ('Azithromycin_AZM', 'Azithromycin'), ('Enrofloxacin_ENR', 'Enrofloxacin'),
    ('Cefoxitin_FOX', 'Cefoxitin'), ('Clindamycin_DA', 'Clindamycin'),
    ('Tetracycline_TE', 'Tetracycline'), ('Gentamicin_CN', 'Gentamicin'),
    ('Chloramphenicol_C', 'Chloramphenicol'), ('TMP_SXT', 'TMP-SXT'),
    ('Nitrofurantoin_F', 'Nitrofurantoin')
]
rows = []
for col, label in agents:
    vc = sa[col].value_counts()
    rows.append((label, int(vc.get('Susceptible', 0)), int(vc.get('Intermediate', 0)), int(vc.get('Resistant', 0))))

ast = pd.DataFrame(rows, columns=['agent', 'Susceptible', 'Intermediate', 'Resistant'])
ast.to_csv(DATA_OUT / 'Figure_S1_AST_profile.csv', index=False)

fig, ax = plt.subplots(figsize=(6.9, 5.1))
y = np.arange(len(ast))
left = np.zeros(len(ast))
for status, color in [('Susceptible', SUS), ('Intermediate', INT), ('Resistant', RES)]:
    vv = 100 * ast[status].to_numpy() / 24
    ax.barh(y, vv, left=left, color=color, edgecolor='white', lw=1, height=0.72, label=status)
    if status == 'Resistant':
        for yy, val, ll in zip(y, vv, left):
            if val >= 10:
                ax.text(ll + val / 2, yy, f'{val:.1f}%', ha='center', va='center', fontsize=7.2, color='white')
    left += vv

ax.set(yticks=y, yticklabels=ast.agent, xlim=(0, 100), xlabel=r'Isolates (%)  |  $n=24$')
ax.set_xticks(range(0, 101, 20), [f'{x}%' for x in range(0, 101, 20)])
ax.invert_yaxis()
ax.grid(axis='x', color=LIGHT, lw=0.55)
ax.set_axisbelow(True)
ax.spines[['top', 'right', 'left']].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.legend(frameon=False, ncol=3, loc='lower center', bbox_to_anchor=(0.5, -0.18))
ax.text(0, 1.035, r'Antimicrobial susceptibility profile of $\it{S.\ aureus}$ isolates', transform=ax.transAxes, fontweight='bold', fontsize=10.5)
fig.tight_layout()
save(fig, SUPP / 'Figure_S1_AST_stacked_profile')

# ==============================================================================
# FIGURE S2: Resistance-class burden by mecA status
# ==============================================================================
print("Generating Figure S2: Resistance-class burden...")
burden.to_csv(DATA_OUT / 'Figure_S2_resistance_class_burden.csv', index=False)
fig, ax = plt.subplots(figsize=(4.7, 3.8))
rng = np.random.default_rng(12)
groups = ['mecA-positive', 'mecA-negative']
cols = [AMBER, BLUE]
for i, (g, c) in enumerate(zip(groups, cols)):
    v = burden[burden.mecA_status.eq(g)].Resistance_classes.to_numpy()
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    x = i
    ax.vlines(x, q1, q3, color=c, lw=5, alpha=0.45, zorder=1)
    ax.hlines(med, x - 0.18, x + 0.18, color='black', lw=1.3, zorder=3)
    ax.scatter(x + rng.uniform(-0.09, 0.09, len(v)), v, s=30, fc=c, ec='black', lw=0.45, zorder=4)
    ax.text(x, 8.48, f'Median {med:.1f}', ha='center', fontsize=8, fontweight='bold', color=c)

ax.axhline(3, color=MID, ls=(0, (3, 2)), lw=0.8)
ax.text(-0.38, 3.12, 'MDR threshold', fontsize=7.2, color=MID)
ax.set(xlim=(-0.48, 1.48), ylim=(0.65, 8.8), xticks=[0, 1], ylabel='Resistant antimicrobial classes')
ax.set_xticklabels([r'$\it{mecA}$ positive' + '\n(n = 18)', r'$\it{mecA}$ negative' + '\n(n = 6)'])
ax.set_yticks(range(1, 9))
ax.grid(axis='y', color=LIGHT, lw=0.55)
ax.set_axisbelow(True)
ax.spines[['top', 'right']].set_visible(False)
ax.text(0, 1.035, r'Resistance-class burden by $\it{mecA}$ status', transform=ax.transAxes, fontweight='bold', fontsize=10.5)
fig.tight_layout()
save(fig, SUPP / 'Figure_S2_resistance_class_burden_by_mecA_status')

# ==============================================================================
# FIGURE S3: Isolate-level AST heatmap
# ==============================================================================
print("Generating Figure S3: Isolate AST heatmap...")
heat_agents = [
    ('Cefoxitin_FOX', 'Cefoxitin'), ('Nitrofurantoin_F', 'Nitrofurantoin'),
    ('TMP_SXT', 'TMP-SXT'), ('Tetracycline_TE', 'Tetracycline'),
    ('Chloramphenicol_C', 'Chloramphenicol'), ('Gentamicin_CN', 'Gentamicin'),
    ('Ciprofloxacin_CIP', 'Ciprofloxacin'), ('Azithromycin_AZM', 'Azithromycin'),
    ('Clindamycin_DA', 'Clindamycin'), ('Erythromycin_E', 'Erythromycin'),
    ('Enrofloxacin_ENR', 'Enrofloxacin')
]
heat = sa.merge(burden[['Sample_ID', 'Resistance_classes']], on='Sample_ID', how='left')

def classify(r):
    if r.PCR_mecA == 'Yes' and r.Cefoxitin_FOX == 'Resistant':
        return 'Typical MRSA'
    if r.PCR_mecA == 'Yes' and r.Cefoxitin_FOX == 'Susceptible':
        return 'OS-MRSA'
    if r.PCR_mecA == 'No' and r.Cefoxitin_FOX == 'Resistant':
        return 'Cefoxitin-resistant\nmecA-negative'
    return 'MSSA'

order = ['Typical MRSA', 'OS-MRSA', 'Cefoxitin-resistant\nmecA-negative', 'MSSA']
heat['classification'] = heat.apply(classify, axis=1)
heat['classification'] = pd.Categorical(heat.classification, categories=order, ordered=True)
heat = heat.sort_values(['classification', 'Resistance_classes', 'Sample_ID'], ascending=[True, False, True]).reset_index(drop=True)
heat[['Sample_ID', 'classification', 'Resistance_classes'] + [x[0] for x in heat_agents]].to_csv(DATA_OUT / 'Figure_S3_isolate_AST_heatmap.csv', index=False)

mapping = {'Susceptible': 0, 'Intermediate': 1, 'Resistant': 2}
arr = np.array([[mapping[heat.loc[i, c]] for c, _ in heat_agents] for i in range(len(heat))])
fig = plt.figure(figsize=(8.1, 7.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 8.5], wspace=0.035)
ag = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1])

ax.imshow(arr, aspect='auto', cmap=ListedColormap([SUS, INT, RES]), norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3))
ax.set_xticks(range(len(heat_agents)), [x[1] for x in heat_agents], rotation=43, ha='right')
ax.set_yticks(range(len(heat)), heat.Sample_ID)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

for i in range(1, len(heat)):
    if heat.loc[i, 'classification'] != heat.loc[i - 1, 'classification']:
        ax.axhline(i - 0.5, color='black', lw=1.1)

classcols = {'Typical MRSA': BLUE, 'OS-MRSA': AMBER, 'Cefoxitin-resistant\nmecA-negative': GREY, 'MSSA': TEAL}
ag.set(xlim=(0, 1), ylim=(len(heat) - 0.5, -0.5), xticks=[], yticks=[])
for s in ag.spines.values():
    s.set_visible(False)

for g in order:
    idx = np.where(heat.classification.astype(str).to_numpy() == g)[0]
    if len(idx):
        ag.add_patch(Rectangle((0.55, idx.min() - 0.5), 0.22, len(idx), fc=classcols[g], ec='white', lw=0.5))
        ag.text(0.48, (idx.min() + idx.max()) / 2, g, ha='right', va='center', fontsize=6.7)

fig.text(0.01, 0.975, 'Isolate-level antimicrobial susceptibility profile', fontweight='bold', fontsize=10.5, va='top')
fig.text(0.01, 0.945, 'Rows grouped by genotype–phenotype classification and ordered by resistance-class burden.', fontsize=7.5, va='top')
ax.legend(handles=[Patch(fc=SUS, label='Susceptible'), Patch(fc=INT, label='Intermediate'), Patch(fc=RES, label='Resistant')], frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.105))
fig.subplots_adjust(left=0.05, right=0.985, top=0.90, bottom=0.17)
save(fig, SUPP / 'Figure_S3_isolate_AST_heatmap')

# ==============================================================================
# FIGURE S4: MecR1 domain architecture schematic
# ==============================================================================
print("Generating Figure S4: MecR1 schematic...")
fig, axes = plt.subplots(2, 1, figsize=(6.8, 3.0), sharex=True, gridspec_kw={'hspace': 0.55})
for ax, title in zip(axes, ['A. Full-length MecR1 reference (585 aa)', 'B. Partial assembly call']):
    ax.set_title(title, loc='left', fontweight='bold', fontsize=9.5, pad=3)
    ax.set(xlim=(0, 600), ylim=(0, 1), yticks=[])
    ax.spines[['left', 'right', 'top']].set_visible(False)

# Reference
axes[0].add_patch(Rectangle((0, 0.27), 328, 0.42, fc='#F6CB56', ec='black', lw=0.65))
axes[0].add_patch(Rectangle((328, 0.27), 257, 0.42, fc='#69A8D1', ec='black', lw=0.65))
axes[0].axvline(204, ymin=0.20, ymax=0.84, color='#A84200', lw=1.6)
axes[0].text(164, 0.48, 'N-terminal transmembrane /\nmetalloprotease region', ha='center', va='center', fontsize=7.8)
axes[0].text(456, 0.48, 'C-terminal β-lactam\nsensor region', ha='center', va='center', fontsize=7.8)
axes[0].text(204, 0.87, 'HExxH 204–208', ha='center', fontsize=7.4, color='#A84200')

# Partial call
axes[1].add_patch(Rectangle((0, 0.27), 328, 0.42, fc='#F6CB56', ec='black', lw=0.65))
axes[1].add_patch(Rectangle((328, 0.27), 257, 0.42, fc='#FAEFEF', ec='#B74545', lw=0.65, hatch='///'))
axes[1].axvline(204, ymin=0.20, ymax=0.84, color='#A84200', lw=1.6)
axes[1].text(164, 0.48, '328-aa predicted ORF\n325 aa aligned to reference', ha='center', va='center', fontsize=7.8)
axes[1].text(456, 0.48, 'Reference-aligned C-terminal\nregion not represented', ha='center', va='center', fontsize=7.8, color='#8E2424')
axes[1].text(204, 0.87, 'HExxH retained', ha='center', fontsize=7.4, color='#A84200')
axes[1].set_xlabel('Amino-acid position in WP_000952923.1 reference')
fig.text(0.5, 0.01, 'Schematic alignment summary; not evidence of the cause of the partial call or cefoxitin phenotype.', ha='center', fontsize=6.7, color=MID)
fig.subplots_adjust(left=0.10, right=0.99, bottom=0.22, top=0.92)
save(fig, SUPP / 'Figure_S4_MecR1_domain_architecture')

print("\nAll publication figures successfully generated in:", OUT)
