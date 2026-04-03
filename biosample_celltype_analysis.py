#!/usr/bin/env python3
"""
Enhanced Cell Type Analysis for Biosamples.csv
==============================================
Extracts base cell type names by stripping donor-specific info (age, sex, treatment),
lists all cell types per organ, and generates focused cell-type plots.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import os
from collections import defaultdict

OUTPUT_DIR = "biosample_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)
sns.set_theme(style="whitegrid", font_scale=1.1)

# ── Load Data ───────────────────────────────────────────────────────────────
df = pd.read_csv("Biosamples.csv")
samples = df[
    df['Biosample'].notna() & (df['Biosample'] != '') &
    df['Organ/Tissue'].notna() & (df['Organ/Tissue'] != '')
].copy()

# ══════════════════════════════════════════════════════════════════════════════
# EXTRACT BASE CELL TYPE NAMES
# ══════════════════════════════════════════════════════════════════════════════
# Many biosample names include donor info like:
#   "CD4-positive, alpha-beta T cell, male adult (38 years) treated with ..."
#   "K562"
#   "GM12878"
# We want to extract the base cell type: "CD4-positive, alpha-beta T cell", "K562", "GM12878"

def extract_base_celltype(name):
    """Strip donor-specific info (sex, age, treatment, disease) from biosample name to get base cell type."""
    if not isinstance(name, str):
        return name

    # Remove treatment info: ", treated with ..."
    name_clean = re.split(r',?\s*treated with\b', name)[0]

    # Remove donor info patterns: ", male/female adult/child/embryo (XX years/days)"
    # Also handles: "(90 or above years)", "(96 days)", "(108 days)"
    # And anything after (disease annotations like "with Alzheimer's disease")
    name_clean = re.split(
        r',?\s*(?:male|female)?\s*(?:adult|child|embryo)\s*\(\d+(?:\s*or\s+above)?\s*(?:years?|days?|weeks?)\).*$',
        name_clean
    )[0]

    # Remove standalone sex/age: ", female", ", male adult"
    name_clean = re.split(r',\s*(?:male|female)(?:\s+adult)?\s*$', name_clean)[0]

    # Remove ", with <disease>" annotations (e.g. ", with multiple sclerosis")
    name_clean = re.split(r',\s*with\s+', name_clean)[0]

    # Remove "nuclear fraction" suffix
    name_clean = re.sub(r',?\s*nuclear fraction\s*$', '', name_clean)

    # Normalize organoid names: "Brain  organoid ( 180 days post differentiation)" -> "Brain organoid"
    name_clean = re.sub(r'\s*\(\s*\d+\s*days?\s+post\s+differentiation\s*\)', '', name_clean)

    # Normalize multiple spaces
    name_clean = re.sub(r'\s{2,}', ' ', name_clean)

    # Remove trailing whitespace and commas
    name_clean = name_clean.strip().rstrip(',').strip()

    return name_clean

samples['base_celltype'] = samples['Biosample'].apply(extract_base_celltype)

# ── Set up dual output: console + text file ──────────────────────────────────
import sys
import io

TEXT_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "celltype_analysis_output.txt")

class TeeOutput:
    """Write to both console and a file simultaneously."""
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.file = open(filepath, 'w', encoding='utf-8')
    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)
    def flush(self):
        self.terminal.flush()
        self.file.flush()
    def close(self):
        self.file.close()

tee = TeeOutput(TEXT_OUTPUT_FILE)
sys.stdout = tee

# ══════════════════════════════════════════════════════════════════════════════
# PRINT COMPREHENSIVE CELL TYPE LISTING
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 80)
print("COMPLETE CELL TYPE LISTING BY ORGAN AND SAMPLE TYPE")
print("=" * 80)

all_celltypes_by_organ = {}

for organ in sorted(samples['Organ/Tissue'].unique()):
    organ_data = samples[samples['Organ/Tissue'] == organ]
    organ_celltypes = {}

    for stype in ['Cell line', 'Primary cell', 'In vitro differentiated cells', 'Organoid', 'Tissue']:
        type_data = organ_data[organ_data['Sample Type'] == stype]
        if len(type_data) == 0:
            continue
        base_types = sorted(type_data['base_celltype'].unique())
        organ_celltypes[stype] = base_types

    all_celltypes_by_organ[organ] = organ_celltypes

# Print the full listing
total_base_celltypes = set()
for organ in sorted(all_celltypes_by_organ.keys()):
    organ_celltypes = all_celltypes_by_organ[organ]
    total_in_organ = sum(len(v) for v in organ_celltypes.values())

    print(f"\n{'─' * 80}")
    print(f"ORGAN: {organ} ({total_in_organ} unique base cell types/samples)")
    print(f"{'─' * 80}")

    for stype in ['Cell line', 'Primary cell', 'In vitro differentiated cells', 'Organoid', 'Tissue']:
        if stype not in organ_celltypes:
            continue
        types = organ_celltypes[stype]
        print(f"\n  [{stype}] ({len(types)} types):")
        for t in types:
            total_base_celltypes.add((organ, stype, t))
            # Count how many individual samples this base type has
            count = len(samples[(samples['Organ/Tissue'] == organ) &
                               (samples['base_celltype'] == t)])
            print(f"    - {t}" + (f"  ({count} samples)" if count > 1 else ""))

print(f"\n\n{'=' * 80}")
print(f"SUMMARY: {len(total_base_celltypes)} total unique (organ, sample_type, base_celltype) combinations")
print(f"{'=' * 80}")

# ── Count base cell types per category ──────────────────────────────────────
all_base = samples.drop_duplicates(['Organ/Tissue', 'Sample Type', 'base_celltype'])
print(f"\nBase cell types by Sample Type:")
for stype in ['Cell line', 'Primary cell', 'In vitro differentiated cells', 'Organoid', 'Tissue']:
    count = len(all_base[all_base['Sample Type'] == stype])
    print(f"  {stype}: {count} unique base types across all organs")

# ── Explicitly list all CELL LINES ──────────────────────────────────────────
print(f"\n\n{'=' * 80}")
print("ALL CELL LINES (base names, grouped by organ)")
print(f"{'=' * 80}")
cl_data = samples[samples['Sample Type'] == 'Cell line']
cl_base = cl_data.drop_duplicates('base_celltype').sort_values(['Organ/Tissue', 'base_celltype'])
cl_by_organ = cl_base.groupby('Organ/Tissue')['base_celltype'].apply(list)
total_cl = 0
for organ, types in cl_by_organ.items():
    print(f"\n  {organ} ({len(types)} cell lines):")
    for t in sorted(types):
        n = len(cl_data[cl_data['base_celltype'] == t])
        print(f"    * {t}" + (f"  [{n} samples]" if n > 1 else ""))
    total_cl += len(types)
print(f"\n  TOTAL unique cell line base names: {total_cl}")

# ── Explicitly list all PRIMARY CELLS ───────────────────────────────────────
print(f"\n\n{'=' * 80}")
print("ALL PRIMARY CELL TYPES (base names, grouped by organ)")
print(f"{'=' * 80}")
pc_data = samples[samples['Sample Type'] == 'Primary cell']
pc_base = pc_data.drop_duplicates(['Organ/Tissue', 'base_celltype']).sort_values(['Organ/Tissue', 'base_celltype'])
pc_by_organ = pc_base.groupby('Organ/Tissue')['base_celltype'].apply(list)
total_pc = 0
for organ, types in pc_by_organ.items():
    print(f"\n  {organ} ({len(types)} primary cell types):")
    for t in sorted(types):
        n = len(pc_data[(pc_data['Organ/Tissue'] == organ) & (pc_data['base_celltype'] == t)])
        print(f"    * {t}" + (f"  [{n} donors/samples]" if n > 1 else ""))
    total_pc += len(types)
print(f"\n  TOTAL unique primary cell base types: {total_pc}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\nGenerating cell-type-focused plots...")

# ── Plot 13: Base Cell Types per Organ (all sample types, stacked) ─────────
base_per_organ = all_base.groupby(['Organ/Tissue', 'Sample Type']).size().unstack(fill_value=0)
base_per_organ['Total'] = base_per_organ.sum(axis=1)
base_per_organ = base_per_organ.sort_values('Total', ascending=True)
base_per_organ = base_per_organ.drop('Total', axis=1)

fig, ax = plt.subplots(figsize=(14, max(10, len(base_per_organ) * 0.4)))
colors = {'Cell line': '#e74c3c', 'Primary cell': '#3498db',
          'In vitro differentiated cells': '#2ecc71', 'Organoid': '#f39c12', 'Tissue': '#95a5a6'}
col_order = [c for c in ['Tissue', 'Primary cell', 'Cell line',
              'In vitro differentiated cells', 'Organoid'] if c in base_per_organ.columns]
base_per_organ[col_order].plot(
    kind='barh', stacked=True, ax=ax,
    color=[colors[c] for c in col_order], edgecolor='white', linewidth=0.5
)
ax.set_xlabel("Number of Unique Base Cell Types", fontsize=13)
ax.set_title("Base Cell Types per Organ/Tissue (by Sample Type)", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type", bbox_to_anchor=(1.02, 1), loc='upper left')
# Add total numbers
for i, (organ, row) in enumerate(base_per_organ.iterrows()):
    total = row.sum()
    ax.text(total + 0.3, i, str(int(total)), va='center', fontsize=9, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/13_base_celltypes_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [13] Base cell types per organ (stacked) saved.")

# ── Plot 14: Cell Lines per Organ (detailed) ──────────────────────────────
cl_base_per_organ = cl_data.groupby('Organ/Tissue')['base_celltype'].nunique().sort_values(ascending=True)
fig, ax = plt.subplots(figsize=(12, max(6, len(cl_base_per_organ) * 0.45)))
bars = ax.barh(range(len(cl_base_per_organ)), cl_base_per_organ.values,
               color=sns.color_palette("Reds_r", len(cl_base_per_organ)))
ax.set_yticks(range(len(cl_base_per_organ)))
ax.set_yticklabels(cl_base_per_organ.index, fontsize=10)
ax.set_xlabel("Number of Unique Cell Lines (base names)", fontsize=13)
ax.set_title("Cell Lines per Organ (Base Names Only)", fontsize=16, fontweight='bold')
for i, v in enumerate(cl_base_per_organ.values):
    ax.text(v + 0.2, i, str(v), va='center', fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/14_cell_lines_base_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [14] Cell lines (base names) per organ saved.")

# ── Plot 15: Primary Cell Types per Organ ─────────────────────────────────
pc_base_per_organ = pc_data.groupby('Organ/Tissue')['base_celltype'].nunique().sort_values(ascending=True)
fig, ax = plt.subplots(figsize=(12, max(6, len(pc_base_per_organ) * 0.45)))
bars = ax.barh(range(len(pc_base_per_organ)), pc_base_per_organ.values,
               color=sns.color_palette("Blues_r", len(pc_base_per_organ)))
ax.set_yticks(range(len(pc_base_per_organ)))
ax.set_yticklabels(pc_base_per_organ.index, fontsize=10)
ax.set_xlabel("Number of Unique Primary Cell Types (base names)", fontsize=13)
ax.set_title("Primary Cell Types per Organ (Base Names Only)", fontsize=16, fontweight='bold')
for i, v in enumerate(pc_base_per_organ.values):
    ax.text(v + 0.2, i, str(v), va='center', fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/15_primary_cells_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [15] Primary cell types per organ saved.")

# ── Plot 16: All Most-Sampled Cell Types (by number of individual samples) ──
base_sample_counts = samples.groupby(['base_celltype', 'Sample Type']).size().reset_index(name='count')
base_totals = base_sample_counts.groupby('base_celltype')['count'].sum().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(14, max(10, len(base_totals) * 0.25)))
pivot = base_sample_counts.pivot_table(index='base_celltype', columns='Sample Type', values='count', fill_value=0)
pivot = pivot.loc[base_totals.index]  # maintain sort order
col_order_present = [c for c in col_order if c in pivot.columns]
pivot[col_order_present].plot(
    kind='barh', stacked=True, ax=ax,
    color=[colors[c] for c in col_order_present], edgecolor='white', linewidth=0.3
)
ax.set_xlabel("Number of Individual Samples", fontsize=13)
ax.set_title("All Cell Types / Tissues by Sample Count", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type", bbox_to_anchor=(1.02, 1), loc='upper left')
ax.invert_yaxis()
for i, (ct, total) in enumerate(base_totals.items()):
    ax.text(total + 0.3, i, str(total), va='center', fontsize=6, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/16_all_celltypes_by_sample_count.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [16] All cell types by sample count saved.")

# ── Plot 17: Blood Cell Types (all) ──────────────────────────────────────────
blood = samples[samples['Organ/Tissue'] == 'Blood']
blood_base = blood.groupby(['base_celltype', 'Sample Type']).size().reset_index(name='count')
blood_totals = blood_base.groupby('base_celltype')['count'].sum().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(14, max(10, len(blood_totals) * 0.3)))
blood_pivot = blood_base.pivot_table(
    index='base_celltype', columns='Sample Type', values='count', fill_value=0
)
blood_pivot = blood_pivot.loc[blood_totals.index]
col_present = [c for c in col_order if c in blood_pivot.columns]
blood_pivot[col_present].plot(
    kind='barh', stacked=True, ax=ax,
    color=[colors[c] for c in col_present], edgecolor='white', linewidth=0.3
)
ax.set_xlabel("Number of Individual Samples", fontsize=13)
ax.set_title("Blood: All Cell Types by Sample Count", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/17_blood_cell_types_all.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [17] Blood all cell types saved.")

# ── Plot 18: Brain Cell Types (all) ──────────────────────────────────────────
brain = samples[samples['Organ/Tissue'] == 'Brain']
brain_base = brain.groupby(['base_celltype', 'Sample Type']).size().reset_index(name='count')
brain_totals = brain_base.groupby('base_celltype')['count'].sum().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(14, max(8, len(brain_totals) * 0.3)))
brain_pivot = brain_base.pivot_table(
    index='base_celltype', columns='Sample Type', values='count', fill_value=0
)
brain_pivot = brain_pivot.loc[brain_totals.index]
col_present = [c for c in col_order if c in brain_pivot.columns]
brain_pivot[col_present].plot(
    kind='barh', stacked=True, ax=ax,
    color=[colors[c] for c in col_present], edgecolor='white', linewidth=0.3
)
ax.set_xlabel("Number of Individual Samples", fontsize=13)
ax.set_title("Brain: All Cell Types by Sample Count", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/18_brain_cell_types_all.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [18] Brain all cell types saved.")

# ── Plot 19: Samples per Donor (how many donors per base cell type) ───────
# For top organs, how many donors contributed to each base cell type?
top6_organs = ['Blood', 'Brain', 'Muscle', 'Kidney', 'Heart', 'Lung']
fig, axes = plt.subplots(2, 3, figsize=(26, 20))
fig.suptitle("Donor/Sample Diversity per Base Cell Type (Top 6 Organs)", fontsize=18, fontweight='bold')

for idx, organ in enumerate(top6_organs):
    ax = axes[idx // 3, idx % 3]
    org_data = samples[samples['Organ/Tissue'] == organ]
    org_counts = org_data.groupby('base_celltype').size().sort_values(ascending=False)

    bar_colors = []
    for ct in org_counts.index:
        stype = org_data[org_data['base_celltype'] == ct]['Sample Type'].iloc[0]
        bar_colors.append(colors.get(stype, '#95a5a6'))

    bars = ax.barh(range(len(org_counts)), org_counts.values, color=bar_colors)
    ax.set_yticks(range(len(org_counts)))
    fontsize = max(4, min(8, 120 // max(len(org_counts), 1)))
    ax.set_yticklabels(org_counts.index, fontsize=fontsize)
    ax.set_xlabel("# Samples")
    ax.set_title(f"{organ} (All {len(org_counts)} types)", fontsize=13, fontweight='bold')
    ax.invert_yaxis()
    for i, v in enumerate(org_counts.values):
        ax.text(v + 0.2, i, str(v), va='center', fontsize=6)

# Add a legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[k], label=k) for k in col_order if k in colors]
fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=11,
           bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig(f"{OUTPUT_DIR}/19_top6_organs_cell_diversity.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [19] Top 6 organs cell diversity saved.")

# ── Plot 20: Organ-CellType Heatmap (top organs x top cell categories) ────
# Categorize base cell types into broader categories for a readable heatmap
def categorize_celltype(name, stype):
    name_lower = name.lower()
    if stype == 'Cell line':
        return 'Cell line'
    if stype == 'Organoid':
        return 'Organoid'
    if stype == 'In vitro differentiated cells':
        return 'In vitro diff.'
    if stype == 'Tissue':
        return 'Tissue'
    # Primary cell subcategories
    if 't cell' in name_lower or 't-cell' in name_lower or 'th1' in name_lower or 'th2' in name_lower or 'th17' in name_lower:
        return 'T cell'
    if 'b cell' in name_lower:
        return 'B cell'
    if 'nk cell' in name_lower or 'natural killer' in name_lower:
        return 'NK cell'
    if 'monocyte' in name_lower:
        return 'Monocyte'
    if 'macrophage' in name_lower:
        return 'Macrophage'
    if 'dendritic' in name_lower:
        return 'Dendritic cell'
    if 'neutrophil' in name_lower:
        return 'Neutrophil'
    if 'fibroblast' in name_lower:
        return 'Fibroblast'
    if 'endothelial' in name_lower:
        return 'Endothelial'
    if 'epithelial' in name_lower or 'keratinocyte' in name_lower:
        return 'Epithelial'
    if 'stem cell' in name_lower or 'progenitor' in name_lower:
        return 'Stem/Progenitor'
    if 'muscle' in name_lower or 'myocyte' in name_lower or 'myoblast' in name_lower:
        return 'Muscle cell'
    if 'neuron' in name_lower or 'astrocyte' in name_lower or 'glial' in name_lower:
        return 'Neural cell'
    return 'Other primary'

samples['cell_category'] = samples.apply(
    lambda row: categorize_celltype(row['base_celltype'], row['Sample Type']), axis=1
)

cat_organ = samples.groupby(['Organ/Tissue', 'cell_category'])['base_celltype'].nunique().unstack(fill_value=0)
# Keep only interesting categories (non-zero in at least a few organs)
cat_cols = [c for c in cat_organ.columns if cat_organ[c].sum() >= 2]
cat_organ = cat_organ[cat_cols]
cat_organ = cat_organ.loc[cat_organ.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(16, max(10, len(cat_organ) * 0.35)))
sns.heatmap(cat_organ, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
            linewidths=0.5, linecolor='white', cbar_kws={'label': 'Unique Base Cell Types'})
ax.set_title("Cell Type Categories per Organ (All Organs)", fontsize=16, fontweight='bold')
ax.set_xlabel("Cell Type Category", fontsize=13)
ax.set_ylabel("")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/20_celltype_category_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [20] Cell type category heatmap saved.")

# ── Plot 21: Treemap-style summary of all cell types ──────────────────────
# Waffle/grid showing organ sizes proportional to cell type count
fig, ax = plt.subplots(figsize=(16, 10))
organ_total = all_base.groupby('Organ/Tissue').size().sort_values(ascending=False)
# Use a grouped horizontal bar showing Tissue vs non-Tissue (cell types)
tissue_counts = all_base[all_base['Sample Type'] == 'Tissue'].groupby('Organ/Tissue').size()
nontissue_counts = all_base[all_base['Sample Type'] != 'Tissue'].groupby('Organ/Tissue').size()

compare_df = pd.DataFrame({
    'Tissue samples': tissue_counts,
    'Cell types (non-tissue)': nontissue_counts
}).fillna(0).astype(int)
compare_df['Total'] = compare_df.sum(axis=1)
compare_df = compare_df.sort_values('Total', ascending=True)
compare_df = compare_df.drop('Total', axis=1)

compare_df.plot(kind='barh', stacked=True, ax=ax,
                color=['#95a5a6', '#e74c3c'], edgecolor='white', linewidth=0.5)
ax.set_xlabel("Number of Unique Base Types", fontsize=13)
ax.set_title("Tissue Samples vs Cell Types per Organ", fontsize=16, fontweight='bold')
ax.legend(title="Category", fontsize=11)
for i, (organ, row) in enumerate(compare_df.iterrows()):
    total = row.sum()
    ax.text(total + 0.3, i, str(int(total)), va='center', fontsize=8, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/21_tissue_vs_celltypes_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [21] Tissue vs cell types per organ saved.")

# ── Final summary stats ───────────────────────────────────────────────────
print(f"\n{'=' * 80}")
print("FINAL CELL TYPE STATISTICS")
print(f"{'=' * 80}")
print(f"  Total biosamples (individual entries):       {len(samples)}")
print(f"  Unique base cell type names (across all):    {samples['base_celltype'].nunique()}")
print(f"  Unique (organ, base_celltype) combos:        {len(samples.drop_duplicates(['Organ/Tissue','base_celltype']))}")
print()
print(f"  By Sample Type (unique base names):")
for stype in ['Cell line', 'Primary cell', 'In vitro differentiated cells', 'Organoid', 'Tissue']:
    n = samples[samples['Sample Type'] == stype]['base_celltype'].nunique()
    print(f"    {stype:<40s} {n:>5d}")
print()
print(f"  By Cell Category:")
for cat in sorted(samples['cell_category'].unique()):
    n = len(samples[samples['cell_category'] == cat])
    n_unique = samples[samples['cell_category'] == cat]['base_celltype'].nunique()
    print(f"    {cat:<30s} {n:>5d} samples, {n_unique:>4d} unique base types")

print(f"\n  Organs with most cell lines:")
for organ, count in cl_base_per_organ.sort_values(ascending=False).head(10).items():
    print(f"    {organ:<35s} {count:>4d} cell lines")

print(f"\n  Organs with most primary cell types:")
for organ, count in pc_base_per_organ.sort_values(ascending=False).head(10).items():
    print(f"    {organ:<35s} {count:>4d} primary cell types")

print(f"\nAll new plots saved to: {OUTPUT_DIR}/")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    if os.path.isfile(fpath):
        fsize = os.path.getsize(fpath) / 1024
        print(f"  {f} ({fsize:.0f} KB)")

# Restore stdout and close the tee file
sys.stdout = tee.terminal
tee.close()
print(f"\nAll output also saved to: {TEXT_OUTPUT_FILE}")
