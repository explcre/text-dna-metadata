#!/usr/bin/env python3
"""
Comprehensive Analysis of ENCODE SCREEN Biosamples.csv
======================================================
Produces statistical summaries and visualizations covering:
 - Cell type / cell line identification and counts
 - Organ/Tissue distribution and biosamples per organ
 - Sample Type breakdown (Tissue, Cell line, Primary cell, etc.)
 - Life Stage distribution
 - Assay coverage analysis (DNase, ATAC, H3K4me3, H3K27ac, CTCF)
 - Collection type distribution
 - RNA-seq availability
 - Assay combinations per biosample
 - Top organs by number of distinct biosamples and cell types
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
import os
import textwrap

# ── Configuration ───────────────────────────────────────────────────────────
OUTPUT_DIR = "biosample_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)
sns.set_theme(style="whitegrid", font_scale=1.1)
PALETTE = sns.color_palette("Set2")

# ── Load Data ───────────────────────────────────────────────────────────────
df = pd.read_csv("Biosamples.csv")
print(f"Total rows in CSV: {len(df)}")
print(f"Columns: {list(df.columns)}\n")

# Filter to actual biosamples (rows that have a Biosample name and an Organ)
samples = df[
    df['Biosample'].notna() & (df['Biosample'] != '') &
    df['Organ/Tissue'].notna() & (df['Organ/Tissue'] != '')
].copy()
print(f"Biosamples with organ annotation: {len(samples)}")
print(f"Unique biosample names: {samples['Biosample'].nunique()}")
print(f"Unique organs/tissues: {samples['Organ/Tissue'].nunique()}\n")

# ── Helper: wrap long labels ────────────────────────────────────────────────
def wrap_labels(labels, width=25):
    return [textwrap.fill(str(l), width) for l in labels]

# ══════════════════════════════════════════════════════════════════════════════
# 1. CELL TYPE / CELL LINE EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. CELL TYPE / CELL LINE ANALYSIS")
print("=" * 70)

cell_lines = samples[samples['Sample Type'] == 'Cell line'].copy()
cell_line_names = sorted(cell_lines['Biosample'].unique())
print(f"\nTotal cell line entries: {len(cell_lines)}")
print(f"Unique cell lines/types: {len(cell_line_names)}")
print("\nAll cell lines/types identified:")
for i, name in enumerate(cell_line_names, 1):
    organ = cell_lines[cell_lines['Biosample'] == name]['Organ/Tissue'].iloc[0]
    print(f"  {i:3d}. {name:<55s} (Organ: {organ})")

# Also identify primary cells
primary_cells = samples[samples['Sample Type'] == 'Primary cell'].copy()
primary_cell_names = sorted(primary_cells['Biosample'].unique())
print(f"\nTotal primary cell entries: {len(primary_cells)}")
print(f"Unique primary cell types: {len(primary_cell_names)}")
print("\nAll primary cell types identified:")
for i, name in enumerate(primary_cell_names, 1):
    organ = primary_cells[primary_cells['Biosample'] == name]['Organ/Tissue'].iloc[0]
    print(f"  {i:3d}. {name:<55s} (Organ: {organ})")

# In-vitro differentiated cells
invitro = samples[samples['Sample Type'] == 'In vitro differentiated cells'].copy()
invitro_names = sorted(invitro['Biosample'].unique())
print(f"\nIn vitro differentiated cell entries: {len(invitro)}")
print(f"Unique in vitro differentiated cell types: {len(invitro_names)}")
for i, name in enumerate(invitro_names, 1):
    organ = invitro[invitro['Biosample'] == name]['Organ/Tissue'].iloc[0]
    print(f"  {i:3d}. {name:<55s} (Organ: {organ})")

# Organoids
organoids = samples[samples['Sample Type'] == 'Organoid'].copy()
organoid_names = sorted(organoids['Biosample'].unique())
print(f"\nOrganoid entries: {len(organoids)}")
print(f"Unique organoid types: {len(organoid_names)}")
for i, name in enumerate(organoid_names, 1):
    organ = organoids[organoids['Biosample'] == name]['Organ/Tissue'].iloc[0]
    print(f"  {i:3d}. {name:<55s} (Organ: {organ})")

# ══════════════════════════════════════════════════════════════════════════════
# 2. PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Plot 1: Sample Type Distribution (pie chart) ───────────────────────────
print("\n\nGenerating plots...")
sample_type_counts = samples['Sample Type'].value_counts()
fig, ax = plt.subplots(figsize=(8, 8))
wedges, texts, autotexts = ax.pie(
    sample_type_counts.values, labels=sample_type_counts.index,
    autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100.*sum(sample_type_counts.values)))})',
    colors=PALETTE, startangle=140, textprops={'fontsize': 12}
)
ax.set_title("Distribution of Sample Types", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/01_sample_type_pie.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [1/12] Sample type pie chart saved.")

# ── Plot 2: Number of Biosamples per Organ (horizontal bar) ────────────────
organ_counts = samples.groupby('Organ/Tissue')['Biosample'].nunique().sort_values(ascending=True)
fig, ax = plt.subplots(figsize=(12, max(10, len(organ_counts) * 0.35)))
bars = ax.barh(range(len(organ_counts)), organ_counts.values, color=sns.color_palette("viridis", len(organ_counts)))
ax.set_yticks(range(len(organ_counts)))
ax.set_yticklabels(organ_counts.index, fontsize=10)
ax.set_xlabel("Number of Unique Biosamples", fontsize=13)
ax.set_title("Number of Unique Biosamples per Organ/Tissue", fontsize=16, fontweight='bold')
for i, v in enumerate(organ_counts.values):
    ax.text(v + 0.3, i, str(v), va='center', fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/02_biosamples_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [2/12] Biosamples per organ chart saved.")

# ── Plot 3: Cell Lines per Organ ────────────────────────────────────────────
cl_per_organ = cell_lines.groupby('Organ/Tissue')['Biosample'].nunique().sort_values(ascending=True)
if len(cl_per_organ) > 0:
    fig, ax = plt.subplots(figsize=(10, max(6, len(cl_per_organ) * 0.4)))
    ax.barh(range(len(cl_per_organ)), cl_per_organ.values, color=sns.color_palette("magma", len(cl_per_organ)))
    ax.set_yticks(range(len(cl_per_organ)))
    ax.set_yticklabels(cl_per_organ.index, fontsize=10)
    ax.set_xlabel("Number of Unique Cell Lines", fontsize=13)
    ax.set_title("Cell Lines per Organ/Tissue", fontsize=16, fontweight='bold')
    for i, v in enumerate(cl_per_organ.values):
        ax.text(v + 0.1, i, str(v), va='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/03_cell_lines_per_organ.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  [3/12] Cell lines per organ chart saved.")

# ── Plot 4: Sample Type Breakdown per Organ (stacked bar) ──────────────────
type_organ = samples.groupby(['Organ/Tissue', 'Sample Type'])['Biosample'].nunique().unstack(fill_value=0)
type_organ = type_organ.loc[type_organ.sum(axis=1).sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(16, 8))
type_organ.plot(kind='bar', stacked=True, ax=ax, colormap='Set2', edgecolor='white', linewidth=0.5)
ax.set_xlabel("Organ/Tissue", fontsize=13)
ax.set_ylabel("Number of Unique Biosamples", fontsize=13)
ax.set_title("Sample Type Composition per Organ/Tissue", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type", bbox_to_anchor=(1.02, 1), loc='upper left')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/04_sample_type_per_organ_stacked.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [4/12] Sample type per organ stacked bar saved.")

# ── Plot 5: Life Stage Distribution ────────────────────────────────────────
life_stage_counts = samples['Life Stage'].value_counts()
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
# Pie
axes[0].pie(life_stage_counts.values, labels=life_stage_counts.index,
            autopct='%1.1f%%', colors=PALETTE[:len(life_stage_counts)], startangle=90)
axes[0].set_title("Life Stage Distribution", fontsize=14, fontweight='bold')
# Life Stage per Organ
ls_organ = samples.groupby(['Organ/Tissue', 'Life Stage'])['Biosample'].nunique().unstack(fill_value=0)
ls_organ = ls_organ.loc[ls_organ.sum(axis=1).sort_values(ascending=False).head(20).index]
ls_organ.plot(kind='barh', stacked=True, ax=axes[1], colormap='coolwarm')
axes[1].set_xlabel("Number of Unique Biosamples")
axes[1].set_title("Life Stage per Organ (Top 20)", fontsize=14, fontweight='bold')
axes[1].legend(title="Life Stage")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/05_life_stage_distribution.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [5/12] Life stage distribution saved.")

# ── Plot 6: Assay Coverage Analysis ────────────────────────────────────────
assay_cols = {
    'DNase': 'DNase Exp. ID',
    'ATAC': 'ATAC Exp. ID',
    'H3K4me3': 'H3K4me3 Exp. ID',
    'H3K27ac': 'H3K27ac Exp. ID',
    'CTCF': 'CTCF Exp. ID'
}

# Count how many biosamples have each assay
assay_presence = {}
for assay_name, col in assay_cols.items():
    has_assay = samples[col].notna() & (samples[col] != '')
    assay_presence[assay_name] = has_assay.sum()

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
# Bar chart of assay counts
assay_names = list(assay_presence.keys())
assay_values = list(assay_presence.values())
colors_assay = sns.color_palette("husl", len(assay_names))
bars = axes[0].bar(assay_names, assay_values, color=colors_assay, edgecolor='black', linewidth=0.5)
axes[0].set_ylabel("Number of Biosamples", fontsize=13)
axes[0].set_title("Assay Coverage Across Biosamples", fontsize=14, fontweight='bold')
for bar, val in zip(bars, assay_values):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                 str(val), ha='center', fontsize=11, fontweight='bold')

# Number of assays per biosample
samples_copy = samples.copy()
samples_copy['num_assays'] = 0
for assay_name, col in assay_cols.items():
    samples_copy['num_assays'] += (samples_copy[col].notna() & (samples_copy[col] != '')).astype(int)
assay_count_dist = samples_copy['num_assays'].value_counts().sort_index()
axes[1].bar(assay_count_dist.index, assay_count_dist.values, color=sns.color_palette("YlOrRd", len(assay_count_dist)),
            edgecolor='black', linewidth=0.5)
axes[1].set_xlabel("Number of Assays per Biosample", fontsize=13)
axes[1].set_ylabel("Count of Biosamples", fontsize=13)
axes[1].set_title("Distribution of Assay Count per Biosample", fontsize=14, fontweight='bold')
axes[1].set_xticks(assay_count_dist.index)
for x, y in zip(assay_count_dist.index, assay_count_dist.values):
    axes[1].text(x, y + 1, str(y), ha='center', fontsize=11)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/06_assay_coverage.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [6/12] Assay coverage analysis saved.")

# ── Plot 7: Assay Availability Heatmap per Organ ──────────────────────────
assay_by_organ = pd.DataFrame()
for assay_name, col in assay_cols.items():
    has = samples[samples[col].notna() & (samples[col] != '')].groupby('Organ/Tissue').size()
    assay_by_organ[assay_name] = has
assay_by_organ = assay_by_organ.fillna(0).astype(int)
assay_by_organ = assay_by_organ.loc[assay_by_organ.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(10, max(10, len(assay_by_organ) * 0.35)))
sns.heatmap(assay_by_organ, annot=True, fmt='d', cmap='YlGnBu', ax=ax,
            linewidths=0.5, linecolor='white', cbar_kws={'label': 'Number of Biosamples'})
ax.set_title("Assay Availability per Organ/Tissue", fontsize=16, fontweight='bold')
ax.set_xlabel("Assay Type", fontsize=13)
ax.set_ylabel("")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/07_assay_heatmap_per_organ.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [7/12] Assay heatmap per organ saved.")

# ── Plot 8: Collection Type Distribution ───────────────────────────────────
collection_counts = samples['Collection'].value_counts()
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
axes[0].pie(collection_counts.values, labels=collection_counts.index,
            autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100.*sum(collection_counts.values)))})',
            colors=PALETTE[:len(collection_counts)], startangle=90)
axes[0].set_title("Collection Type Distribution", fontsize=14, fontweight='bold')

# Collection per organ (top 20)
coll_organ = samples.groupby(['Organ/Tissue', 'Collection'])['Biosample'].nunique().unstack(fill_value=0)
coll_organ = coll_organ.loc[coll_organ.sum(axis=1).sort_values(ascending=False).head(20).index]
coll_organ.plot(kind='barh', stacked=True, ax=axes[1], colormap='Pastel1')
axes[1].set_xlabel("Number of Unique Biosamples")
axes[1].set_title("Collection Type per Organ (Top 20)", fontsize=14, fontweight='bold')
axes[1].legend(title="Collection")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/08_collection_distribution.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [8/12] Collection distribution saved.")

# ── Plot 9: RNA-seq Availability ───────────────────────────────────────────
rna_counts = samples['Has RNA Seq'].value_counts()
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
axes[0].pie(rna_counts.values, labels=['Has RNA-seq' if x == 'yes' else 'No RNA-seq' for x in rna_counts.index],
            autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100.*sum(rna_counts.values)))})',
            colors=['#2ecc71', '#e74c3c'], startangle=90, explode=[0.05]*len(rna_counts))
axes[0].set_title("RNA-seq Availability", fontsize=14, fontweight='bold')

# RNA-seq per organ
rna_organ = samples.groupby(['Organ/Tissue', 'Has RNA Seq'])['Biosample'].nunique().unstack(fill_value=0)
if 'yes' in rna_organ.columns:
    rna_organ = rna_organ.sort_values('yes', ascending=True)
rna_organ.plot(kind='barh', stacked=True, ax=axes[1], color=['#e74c3c', '#2ecc71'])
axes[1].set_xlabel("Number of Unique Biosamples")
axes[1].set_title("RNA-seq Availability per Organ", fontsize=14, fontweight='bold')
axes[1].legend(['No RNA-seq', 'Has RNA-seq'])
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/09_rna_seq_availability.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [9/12] RNA-seq availability saved.")

# ── Plot 10: Top 15 Organs with Most Cell Types (all sample types) ─────────
# Show unique biosample count per organ, colored by sample type
top_organs = organ_counts.sort_values(ascending=False).head(15).index.tolist()
top_organ_data = type_organ.loc[top_organs]

fig, ax = plt.subplots(figsize=(14, 8))
top_organ_data.plot(kind='bar', stacked=True, ax=ax, colormap='tab10', edgecolor='white', linewidth=0.5)
ax.set_ylabel("Number of Unique Biosamples", fontsize=13)
ax.set_xlabel("")
ax.set_title("Top 15 Organs by Biosample Diversity (by Sample Type)", fontsize=16, fontweight='bold')
ax.legend(title="Sample Type", bbox_to_anchor=(1.02, 1), loc='upper left')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
# Add total count on top
for i, organ in enumerate(top_organs):
    total = type_organ.loc[organ].sum()
    ax.text(i, total + 0.5, str(total), ha='center', fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/10_top15_organs_diversity.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [10/12] Top 15 organs diversity saved.")

# ── Plot 11: Assay Combination Analysis ────────────────────────────────────
# What combinations of assays are most common?
samples_copy['assay_combo'] = ''
for assay_name, col in assay_cols.items():
    mask = samples_copy[col].notna() & (samples_copy[col] != '')
    samples_copy.loc[mask, 'assay_combo'] = samples_copy.loc[mask, 'assay_combo'] + assay_name + '+'
samples_copy['assay_combo'] = samples_copy['assay_combo'].str.rstrip('+')
samples_copy.loc[samples_copy['assay_combo'] == '', 'assay_combo'] = 'None'

combo_counts = samples_copy['assay_combo'].value_counts().head(15)
fig, ax = plt.subplots(figsize=(12, 7))
bars = ax.barh(range(len(combo_counts)), combo_counts.values,
               color=sns.color_palette("Spectral", len(combo_counts)))
ax.set_yticks(range(len(combo_counts)))
ax.set_yticklabels(combo_counts.index, fontsize=10)
ax.set_xlabel("Number of Biosamples", fontsize=13)
ax.set_title("Top 15 Assay Combinations Across Biosamples", fontsize=16, fontweight='bold')
for i, v in enumerate(combo_counts.values):
    ax.text(v + 0.5, i, str(v), va='center', fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/11_assay_combinations.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [11/12] Assay combinations saved.")

# ── Plot 12: Comprehensive Summary Dashboard ──────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(20, 14))
fig.suptitle("ENCODE SCREEN Biosamples - Comprehensive Dashboard", fontsize=20, fontweight='bold', y=1.01)

# (a) Sample Type
axes[0, 0].pie(sample_type_counts.values, labels=sample_type_counts.index,
               autopct='%1.1f%%', colors=PALETTE, startangle=140)
axes[0, 0].set_title("Sample Type Distribution", fontweight='bold')

# (b) Top 10 Organs
top10 = organ_counts.sort_values(ascending=False).head(10)
axes[0, 1].barh(range(len(top10)), top10.values, color=sns.color_palette("viridis", len(top10)))
axes[0, 1].set_yticks(range(len(top10)))
axes[0, 1].set_yticklabels(top10.index, fontsize=9)
axes[0, 1].set_xlabel("Unique Biosamples")
axes[0, 1].set_title("Top 10 Organs by Biosample Count", fontweight='bold')
axes[0, 1].invert_yaxis()

# (c) Assay Coverage
axes[0, 2].bar(assay_names, assay_values, color=colors_assay, edgecolor='black', linewidth=0.5)
axes[0, 2].set_ylabel("Biosamples")
axes[0, 2].set_title("Assay Coverage", fontweight='bold')

# (d) Life Stage
axes[1, 0].pie(life_stage_counts.values, labels=life_stage_counts.index,
               autopct='%1.1f%%', colors=['#3498db', '#e67e22'], startangle=90)
axes[1, 0].set_title("Life Stage Distribution", fontweight='bold')

# (e) Collection Type
axes[1, 1].pie(collection_counts.values, labels=collection_counts.index,
               autopct='%1.1f%%', colors=PALETTE[:len(collection_counts)], startangle=90)
axes[1, 1].set_title("Collection Type", fontweight='bold')

# (f) RNA-seq
axes[1, 2].pie(rna_counts.values,
               labels=['Has RNA-seq' if x == 'yes' else 'No RNA-seq' for x in rna_counts.index],
               autopct='%1.1f%%', colors=['#2ecc71', '#e74c3c'], startangle=90)
axes[1, 2].set_title("RNA-seq Availability", fontweight='bold')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/12_comprehensive_dashboard.png", dpi=150, bbox_inches='tight')
plt.close()
print("  [12/12] Comprehensive dashboard saved.")

# ══════════════════════════════════════════════════════════════════════════════
# 3. STATISTICAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STATISTICAL SUMMARY")
print("=" * 70)

print(f"\n{'Metric':<50s} {'Value':>10s}")
print("-" * 62)
print(f"{'Total biosample entries':<50s} {len(samples):>10d}")
print(f"{'Unique biosample names':<50s} {samples['Biosample'].nunique():>10d}")
print(f"{'Unique organs/tissues':<50s} {samples['Organ/Tissue'].nunique():>10d}")
print(f"{'Unique cell lines':<50s} {len(cell_line_names):>10d}")
print(f"{'Unique primary cell types':<50s} {len(primary_cell_names):>10d}")
print(f"{'Unique in vitro differentiated cells':<50s} {len(invitro_names):>10d}")
print(f"{'Unique organoid types':<50s} {len(organoid_names):>10d}")
tissue_count = samples[samples['Sample Type'] == 'Tissue']['Biosample'].nunique()
print(f"{'Unique tissue samples':<50s} {tissue_count:>10d}")

print(f"\n--- Assay Coverage ---")
for assay_name, count in assay_presence.items():
    pct = count / len(samples) * 100
    print(f"{'Biosamples with ' + assay_name:<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- Assays per Biosample ---")
for n_assays in sorted(assay_count_dist.index):
    count = assay_count_dist[n_assays]
    pct = count / len(samples) * 100
    print(f"{'Biosamples with ' + str(n_assays) + ' assay(s)':<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- RNA-seq ---")
for val, count in rna_counts.items():
    pct = count / len(samples) * 100
    label = "With RNA-seq" if val == 'yes' else "Without RNA-seq"
    print(f"{label:<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- Sample Type ---")
for stype, count in sample_type_counts.items():
    pct = count / len(samples) * 100
    print(f"{stype:<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- Life Stage ---")
for ls, count in life_stage_counts.items():
    pct = count / len(samples) * 100
    print(f"{ls:<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- Collection Type ---")
for coll, count in collection_counts.items():
    pct = count / len(samples) * 100
    print(f"{coll:<50s} {count:>6d} ({pct:.1f}%)")

print(f"\n--- Organs with Most Biosamples (Top 10) ---")
top10_organs = organ_counts.sort_values(ascending=False).head(10)
for organ, count in top10_organs.items():
    print(f"  {organ:<45s} {count:>5d} biosamples")

print(f"\n--- Organs with Most Cell Lines (Top 10) ---")
cl_top = cl_per_organ.sort_values(ascending=False).head(10)
for organ, count in cl_top.items():
    print(f"  {organ:<45s} {count:>5d} cell lines")

print(f"\n\nAll plots saved to: {OUTPUT_DIR}/")
print("Files generated:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fsize = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / 1024
    print(f"  {f} ({fsize:.0f} KB)")
