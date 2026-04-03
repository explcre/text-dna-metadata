#!/usr/bin/env python3
"""
Comprehensive analysis of cell type mapping between:
  - entex-files.txt (ENTEx download URLs)
  - entex-metadata.tsv (ENTEx experiment metadata)
  - Biosamples.csv (SCREEN biosample registry)

Outputs: analysis text report + plots
"""

import os
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
from difflib import SequenceMatcher

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "celltype_mapping_output")
os.makedirs(OUT_DIR, exist_ok=True)

REPORT_FILE = os.path.join(OUT_DIR, "celltype_mapping_report.txt")
report_lines = []

def report(text=""):
    report_lines.append(text)
    print(text)

def save_report():
    with open(REPORT_FILE, 'w') as f:
        f.write('\n'.join(report_lines))
    print(f"\n[Report saved to {REPORT_FILE}]")

def savefig(name):
    path = os.path.join(OUT_DIR, name)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    report(f"  [Plot saved: {name}]")

# ============================================================
# 1. LOAD DATA
# ============================================================
BASE = os.path.dirname(os.path.abspath(__file__))

report("=" * 80)
report("COMPREHENSIVE CELL TYPE MAPPING ANALYSIS")
report("ENTEx files / ENTEx metadata / Biosamples.csv")
report("=" * 80)

# --- entex-files.txt ---
with open(os.path.join(BASE, "entex-files.txt")) as f:
    raw_lines = [l.strip().strip('"') for l in f if l.strip()]

entex_urls = [l for l in raw_lines if l.startswith("https://www.encodeproject.org/files/")]
entex_file_ids = []
for url in entex_urls:
    m = re.search(r'/files/(ENCFF\w+)/', url)
    if m:
        entex_file_ids.append(m.group(1))
entex_file_ids_set = set(entex_file_ids)

report(f"\n--- entex-files.txt ---")
report(f"Total lines: {len(raw_lines)}")
report(f"File download URLs: {len(entex_urls)}")
report(f"Unique ENCFF IDs extracted: {len(entex_file_ids_set)}")

# file extensions
ext_counter = Counter()
for url in entex_urls:
    if '.bigWig' in url:
        ext_counter['bigWig'] += 1
    elif '.bed.gz' in url:
        ext_counter['bed.gz'] += 1
    elif '.bigBed' in url:
        ext_counter['bigBed'] += 1
    elif '.bam' in url:
        ext_counter['bam'] += 1
    else:
        ext_counter['other'] += 1
report(f"File types in entex-files.txt: {dict(ext_counter)}")

# --- entex-metadata.tsv ---
meta = pd.read_csv(os.path.join(BASE, "entex-metadata.tsv"), sep='\t', low_memory=False)
report(f"\n--- entex-metadata.tsv ---")
report(f"Rows: {len(meta)}, Columns: {len(meta.columns)}")
report(f"Unique File accessions: {meta['File accession'].nunique()}")
report(f"Unique Experiment accessions: {meta['Experiment accession'].nunique()}")
report(f"Unique Biosample term names: {meta['Biosample term name'].nunique()}")
report(f"Biosample term names: {sorted(meta['Biosample term name'].dropna().unique())}")
report(f"Unique Biosample term ids: {meta['Biosample term id'].nunique()}")
report(f"Assays: {sorted(meta['Assay'].dropna().unique())}")
report(f"Donors: {sorted(meta['Donor(s)'].dropna().unique())}")

# --- Biosamples.csv ---
bio = pd.read_csv(os.path.join(BASE, "Biosamples.csv"), low_memory=False)
report(f"\n--- Biosamples.csv ---")
report(f"Rows: {len(bio)}, Columns: {len(bio.columns)}")
report(f"Unique Tissue/Biosample: {bio['Tissue/Biosample'].nunique()}")
report(f"Unique Biosample: {bio['Biosample'].nunique()}")
report(f"Unique Organ/Tissue: {bio['Organ/Tissue'].nunique()}")
report(f"Organ/Tissue values: {sorted(bio['Organ/Tissue'].dropna().unique())}")

# ============================================================
# 2. MAPPING ANALYSIS: ENTEX-FILES ↔ ENTEX-METADATA (File accession)
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 1: entex-files.txt ENCFF IDs ↔ entex-metadata File accessions")
report("=" * 80)

meta_file_ids = set(meta['File accession'].dropna().unique())
overlap_file = entex_file_ids_set & meta_file_ids
only_entex = entex_file_ids_set - meta_file_ids
only_meta = meta_file_ids - entex_file_ids_set

report(f"ENCFF IDs in entex-files.txt: {len(entex_file_ids_set)}")
report(f"File accessions in entex-metadata: {len(meta_file_ids)}")
report(f"Overlap (matched): {len(overlap_file)}")
report(f"Only in entex-files.txt: {len(only_entex)}")
report(f"Only in entex-metadata: {len(only_meta)}")
report(f"Match rate (entex-files → metadata): {len(overlap_file)/len(entex_file_ids_set)*100:.1f}%")

# ============================================================
# 3. MAPPING ANALYSIS: BIOSAMPLE TERM NAME ↔ TISSUE/BIOSAMPLE
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 2: Biosample term name (entex-metadata) ↔ Tissue/Biosample (Biosamples.csv)")
report("=" * 80)

entex_terms = sorted(meta['Biosample term name'].dropna().unique())
bio_tissues = sorted(bio['Tissue/Biosample'].dropna().unique())
bio_biosamples = sorted(bio['Biosample'].dropna().unique())
bio_organs = sorted(bio['Organ/Tissue'].dropna().unique())

# Exact matches
exact_tissue = set(entex_terms) & set(bio_tissues)
exact_biosample = set(entex_terms) & set(bio_biosamples)
exact_organ = set(entex_terms) & set(bio_organs)

report(f"\nENTEx unique biosample term names: {len(entex_terms)}")
report(f"Biosamples.csv unique Tissue/Biosample: {len(bio_tissues)}")
report(f"Biosamples.csv unique Biosample: {len(bio_biosamples)}")
report(f"Biosamples.csv unique Organ/Tissue: {len(bio_organs)}")

report(f"\nExact matches (term name ↔ Tissue/Biosample): {len(exact_tissue)}")
if exact_tissue:
    report(f"  Matched: {sorted(exact_tissue)}")

report(f"Exact matches (term name ↔ Biosample): {len(exact_biosample)}")
if exact_biosample:
    report(f"  Matched: {sorted(exact_biosample)}")

report(f"Exact matches (term name ↔ Organ/Tissue): {len(exact_organ)}")
if exact_organ:
    report(f"  Matched: {sorted(exact_organ)}")

# Case-insensitive / substring matching
report(f"\n--- Fuzzy / Substring Matching (term name → Tissue/Biosample) ---")
entex_lower = {t: t.lower() for t in entex_terms}
bio_tissue_lower = {t: t.lower() for t in bio_tissues}

fuzzy_matches = defaultdict(list)
for et, et_low in entex_lower.items():
    for bt, bt_low in bio_tissue_lower.items():
        # case-insensitive exact
        if et_low == bt_low:
            fuzzy_matches[et].append(('exact_ci', bt))
        # entex term is substring of biosample tissue
        elif et_low in bt_low:
            fuzzy_matches[et].append(('entex_substring_of_bio', bt))
        # biosample tissue is substring of entex term
        elif bt_low in et_low and len(bt_low) > 3:
            fuzzy_matches[et].append(('bio_substring_of_entex', bt))

report(f"ENTEx terms with at least one fuzzy match: {len(fuzzy_matches)}/{len(entex_terms)}")
for et in sorted(fuzzy_matches.keys()):
    matches = fuzzy_matches[et]
    report(f"  '{et}' →")
    for mtype, bt in matches[:5]:
        report(f"    [{mtype}] '{bt}'")
    if len(matches) > 5:
        report(f"    ... and {len(matches)-5} more")

# Also try matching against Organ/Tissue
report(f"\n--- Fuzzy / Substring Matching (term name → Organ/Tissue) ---")
bio_organ_lower = {t: t.lower() for t in bio_organs}
organ_matches = defaultdict(list)
for et, et_low in entex_lower.items():
    for bo, bo_low in bio_organ_lower.items():
        if et_low == bo_low:
            organ_matches[et].append(('exact_ci', bo))
        elif et_low in bo_low:
            organ_matches[et].append(('entex_substring', bo))
        elif bo_low in et_low and len(bo_low) > 3:
            organ_matches[et].append(('organ_substring', bo))

report(f"ENTEx terms with Organ/Tissue match: {len(organ_matches)}/{len(entex_terms)}")
for et in sorted(organ_matches.keys()):
    matches = organ_matches[et]
    report(f"  '{et}' → {[(m, b) for m, b in matches[:5]]}")

# Sequence similarity top matches
report(f"\n--- Top Sequence Similarity Matches (term name → Biosample) ---")
similarity_results = []
for et in entex_terms:
    best_score = 0
    best_match = ""
    for bt in bio_tissues:
        score = SequenceMatcher(None, et.lower(), bt.lower()).ratio()
        if score > best_score:
            best_score = score
            best_match = bt
    similarity_results.append((et, best_match, best_score))
    report(f"  '{et}' → '{best_match}' (similarity: {best_score:.3f})")

# ============================================================
# 4. MAPPING: EXPERIMENT IDs (ENCSR*)
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 3: Experiment accession IDs (ENCSR*)")
report("=" * 80)

meta_exp_ids = set(meta['Experiment accession'].dropna().unique())

# Extract all ENCSR IDs from Biosamples.csv
bio_exp_cols = ['DNase Exp. ID', 'ATAC Exp. ID', 'H3K4me3 Exp. ID', 'H3K27ac Exp. ID', 'CTCF Exp. ID']
bio_exp_ids = set()
for col in bio_exp_cols:
    if col in bio.columns:
        bio_exp_ids.update(bio[col].dropna().unique())

overlap_exp = meta_exp_ids & bio_exp_ids
report(f"Experiment IDs in entex-metadata: {len(meta_exp_ids)}")
report(f"Experiment IDs in Biosamples.csv: {len(bio_exp_ids)}")
report(f"Overlapping Experiment IDs: {len(overlap_exp)}")
report(f"Match rate (entex → bio): {len(overlap_exp)/max(len(meta_exp_ids),1)*100:.1f}%")
report(f"Match rate (bio → entex): {len(overlap_exp)/max(len(bio_exp_ids),1)*100:.1f}%")

if overlap_exp:
    report(f"\nSample overlapping Experiment IDs (first 20): {sorted(overlap_exp)[:20]}")
    # For overlapping experiment IDs, show what cell types they connect
    report(f"\n--- Cell type connections via shared Experiment IDs ---")
    count = 0
    for eid in sorted(overlap_exp)[:30]:
        meta_rows = meta[meta['Experiment accession'] == eid]
        entex_term = meta_rows['Biosample term name'].iloc[0] if len(meta_rows) > 0 else 'N/A'
        entex_assay = meta_rows['Assay'].iloc[0] if len(meta_rows) > 0 else 'N/A'
        bio_matches = []
        for col in bio_exp_cols:
            if col in bio.columns:
                matching = bio[bio[col] == eid]
                if len(matching) > 0:
                    for _, row in matching.iterrows():
                        bio_matches.append((col, row.get('Tissue/Biosample', 'N/A'), row.get('Organ/Tissue', 'N/A')))
        if bio_matches:
            report(f"  {eid}: ENTEx='{entex_term}' ({entex_assay})")
            for col, tissue, organ in bio_matches:
                report(f"    → Biosamples: '{tissue}' (Organ: '{organ}') via {col}")
            count += 1
    report(f"  (Showed {count} connections)")

# ============================================================
# 5. MAPPING: FILE IDs (ENCFF*) across all tables
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 4: File accession IDs (ENCFF*) across all three sources")
report("=" * 80)

# Extract ENCFF from Biosamples.csv
bio_file_cols = ['DNase File ID', 'ATAC File ID', 'H3K4me3 File ID', 'H3K27ac File ID', 'CTCF File ID']
bio_file_ids = set()
for col in bio_file_cols:
    if col in bio.columns:
        bio_file_ids.update(bio[col].dropna().unique())

# Also extract ENCFF IDs from cCRE bed/bigBed URLs
bio_ccre_file_ids = set()
for col in ['cCREs (.bed)', 'cCREs (.bigBed)']:
    if col in bio.columns:
        for val in bio[col].dropna():
            matches = re.findall(r'(ENCFF\w+)', str(val))
            bio_ccre_file_ids.update(matches)

report(f"ENCFF IDs in entex-files.txt: {len(entex_file_ids_set)}")
report(f"ENCFF IDs in entex-metadata (File accession): {len(meta_file_ids)}")
report(f"ENCFF IDs in Biosamples.csv (assay file columns): {len(bio_file_ids)}")
report(f"ENCFF IDs in Biosamples.csv (cCRE URLs): {len(bio_ccre_file_ids)}")

all_bio_file = bio_file_ids | bio_ccre_file_ids
overlap_entex_bio = entex_file_ids_set & all_bio_file
overlap_meta_bio = meta_file_ids & all_bio_file
triple_overlap = entex_file_ids_set & meta_file_ids & all_bio_file

report(f"\nOverlap entex-files ↔ Biosamples.csv: {len(overlap_entex_bio)}")
report(f"Overlap entex-metadata ↔ Biosamples.csv: {len(overlap_meta_bio)}")
report(f"Triple overlap (all three): {len(triple_overlap)}")

if overlap_meta_bio:
    report(f"\nSample ENCFF IDs shared between entex-metadata and Biosamples.csv (first 15):")
    for fid in sorted(overlap_meta_bio)[:15]:
        meta_row = meta[meta['File accession'] == fid].iloc[0] if fid in meta['File accession'].values else None
        if meta_row is not None:
            report(f"  {fid}: ENTEx term='{meta_row['Biosample term name']}', Assay='{meta_row['Assay']}'")
        # find in bio
        for col in bio_file_cols:
            if col in bio.columns:
                brows = bio[bio[col] == fid]
                for _, br in brows.iterrows():
                    report(f"    → Biosamples: '{br['Tissue/Biosample']}' (Organ: '{br.get('Organ/Tissue','N/A')}') via {col}")

# ============================================================
# 6. MAPPING: Biosample term id (UBERON/CL ontology)
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 5: Biosample term id (ontology IDs like UBERON:*, CL:*)")
report("=" * 80)

term_ids = meta[['Biosample term id', 'Biosample term name']].drop_duplicates()
report(f"Unique (term_id, term_name) pairs: {len(term_ids)}")
for _, row in term_ids.sort_values('Biosample term name').iterrows():
    report(f"  {row['Biosample term id']} → {row['Biosample term name']}")

# ============================================================
# 7. ASSAY COVERAGE COMPARISON
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 6: Assay type comparison")
report("=" * 80)

entex_assays = sorted(meta['Assay'].dropna().unique())
bio_assay_types = []
for col in bio.columns:
    if 'Exp. ID' in col:
        assay = col.replace(' Exp. ID', '')
        bio_assay_types.append(assay)

report(f"ENTEx assays: {entex_assays}")
report(f"Biosamples.csv assay types (from columns): {bio_assay_types}")

# Per assay-tissue breakdown in ENTEx
report(f"\n--- ENTEx: Assay × Biosample term name matrix ---")
assay_tissue = meta.groupby(['Assay', 'Biosample term name']).size().reset_index(name='count')
pivot = assay_tissue.pivot(index='Biosample term name', columns='Assay', values='count').fillna(0).astype(int)
report(pivot.to_string())

# ============================================================
# 8. DONOR ANALYSIS
# ============================================================
report("\n" + "=" * 80)
report("MAPPING 7: Donor analysis in ENTEx metadata")
report("=" * 80)

donors = meta['Donor(s)'].dropna().unique()
report(f"Unique donors: {len(donors)}")
for d in sorted(donors):
    tissues = meta[meta['Donor(s)'] == d]['Biosample term name'].unique()
    report(f"  {d}: {len(tissues)} tissues → {sorted(tissues)}")

# ============================================================
# 9. COMPREHENSIVE MAPPING SUMMARY TABLE
# ============================================================
report("\n" + "=" * 80)
report("COMPREHENSIVE MAPPING SUMMARY")
report("=" * 80)

report(f"""
Linkage Keys Between Tables:
┌─────────────────────────┬───────────────────────────┬──────────────────────────┐
│ entex-files.txt         │ entex-metadata.tsv        │ Biosamples.csv           │
├─────────────────────────┼───────────────────────────┼──────────────────────────┤
│ ENCFF* in URL           │ File accession (ENCFF*)   │ DNase/ATAC/... File ID   │
│                         │                           │ cCRE URLs contain ENCFF  │
│                         │ Experiment accession      │ DNase/ATAC/... Exp. ID   │
│                         │ (ENCSR*)                  │ (ENCSR*)                 │
│                         │ Biosample term name       │ Tissue/Biosample         │
│                         │ (e.g. 'ovary')            │ (e.g. 'Ovary, ...')      │
│                         │ Biosample term name       │ Organ/Tissue             │
│                         │ (normalized tissue name)  │ (organ category)         │
│                         │ Assay                     │ Column-based assay types │
│                         │ (WGBS, RRBS, ...)         │ (DNase, ATAC, H3K4me3,  │
│                         │                           │  H3K27ac, CTCF)          │
│                         │ File download URL         │ Signal URLs              │
│                         │ (matches entex-files URLs)│ (encodeproject.org URLs) │
└─────────────────────────┴───────────────────────────┴──────────────────────────┘

Match Statistics:
  entex-files ↔ entex-metadata (ENCFF):   {len(overlap_file):>5} matches ({len(overlap_file)/max(len(entex_file_ids_set),1)*100:.1f}%)
  entex-metadata ↔ Biosamples (ENCSR):    {len(overlap_exp):>5} matches
  entex-metadata ↔ Biosamples (ENCFF):    {len(overlap_meta_bio):>5} matches
  entex-files ↔ Biosamples (ENCFF):       {len(overlap_entex_bio):>5} matches
  Triple overlap (all three, ENCFF):      {len(triple_overlap):>5} matches
  Cell type name exact matches:           {len(exact_tissue):>5} (term name ↔ Tissue/Biosample)
  Cell type → Organ/Tissue matches:       {len(organ_matches):>5} fuzzy matches
""")

# ============================================================
# 10. UNMATCHED ANALYSIS
# ============================================================
report("=" * 80)
report("UNMATCHED CELL TYPES ANALYSIS")
report("=" * 80)

matched_entex_terms = set(fuzzy_matches.keys()) | set(organ_matches.keys())
unmatched = set(entex_terms) - matched_entex_terms
report(f"\nENTEx biosample terms with NO match in Biosamples.csv: {len(unmatched)}")
for t in sorted(unmatched):
    report(f"  - {t}")

# Biosamples organs with no match in ENTEx
bio_organs_lower_set = {o.lower() for o in bio_organs}
entex_terms_lower_set = {t.lower() for t in entex_terms}
unmatched_organs = [o for o in bio_organs if o.lower() not in entex_terms_lower_set]
report(f"\nBiosamples.csv Organ/Tissue values with no exact match in ENTEx: {len(unmatched_organs)}")
for o in sorted(unmatched_organs)[:30]:
    report(f"  - {o}")
if len(unmatched_organs) > 30:
    report(f"  ... and {len(unmatched_organs)-30} more")


# ============================================================
# PLOTS
# ============================================================
report("\n" + "=" * 80)
report("GENERATING PLOTS")
report("=" * 80)

sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 9})

# --- Plot 1: File type distribution in entex-files.txt ---
fig, ax = plt.subplots(figsize=(6, 4))
labels, values = zip(*ext_counter.most_common())
ax.bar(labels, values, color=sns.color_palette("Set2"))
ax.set_title("File Types in entex-files.txt")
ax.set_ylabel("Count")
for i, v in enumerate(values):
    ax.text(i, v + 5, str(v), ha='center', fontsize=9)
savefig("01_entex_files_type_distribution.png")

# --- Plot 2: Venn-like overlap of ENCFF IDs ---
fig, ax = plt.subplots(figsize=(8, 5))
categories = ['entex-files\nonly', 'entex-meta\nonly', 'Biosamples\nonly',
              'files∩meta', 'files∩bio', 'meta∩bio', 'All three']
vals = [
    len(entex_file_ids_set - meta_file_ids - all_bio_file),
    len(meta_file_ids - entex_file_ids_set - all_bio_file),
    len(all_bio_file - entex_file_ids_set - meta_file_ids),
    len((entex_file_ids_set & meta_file_ids) - all_bio_file),
    len((entex_file_ids_set & all_bio_file) - meta_file_ids),
    len((meta_file_ids & all_bio_file) - entex_file_ids_set),
    len(triple_overlap)
]
colors = sns.color_palette("Set3", len(categories))
bars = ax.bar(categories, vals, color=colors)
ax.set_title("ENCFF ID Overlap Across Three Data Sources")
ax.set_ylabel("Count")
for bar, v in zip(bars, vals):
    if v > 0:
        ax.text(bar.get_x() + bar.get_width()/2, v + 5, str(v), ha='center', fontsize=8)
savefig("02_encff_id_overlap.png")

# --- Plot 3: Experiment ID overlap ---
fig, ax = plt.subplots(figsize=(6, 4))
cats = ['ENTEx only', 'Biosamples only', 'Shared']
vals3 = [len(meta_exp_ids - bio_exp_ids), len(bio_exp_ids - meta_exp_ids), len(overlap_exp)]
ax.bar(cats, vals3, color=['#ff9999', '#99ccff', '#99ff99'])
ax.set_title("Experiment ID (ENCSR*) Overlap")
ax.set_ylabel("Count")
for i, v in enumerate(vals3):
    ax.text(i, v + 2, str(v), ha='center')
savefig("03_experiment_id_overlap.png")

# --- Plot 4: ENTEx assay distribution ---
fig, ax = plt.subplots(figsize=(8, 5))
assay_counts = meta['Assay'].value_counts()
assay_counts.plot(kind='barh', ax=ax, color=sns.color_palette("husl", len(assay_counts)))
ax.set_title("ENTEx Assay Distribution (entex-metadata.tsv)")
ax.set_xlabel("Number of Files")
savefig("04_entex_assay_distribution.png")

# --- Plot 5: Biosample term name distribution in ENTEx ---
fig, ax = plt.subplots(figsize=(10, 8))
term_counts = meta['Biosample term name'].value_counts()
term_counts.plot(kind='barh', ax=ax, color=sns.color_palette("tab20", len(term_counts)))
ax.set_title("ENTEx Biosample Term Name Distribution")
ax.set_xlabel("Number of Files")
savefig("05_entex_biosample_term_distribution.png")

# --- Plot 6: Assay × Tissue heatmap ---
fig, ax = plt.subplots(figsize=(12, 10))
pivot_plot = pivot.copy()
sns.heatmap(pivot_plot, annot=True, fmt='d', cmap='YlOrRd', ax=ax, linewidths=0.5)
ax.set_title("ENTEx: File Count by Assay × Biosample Term Name")
savefig("06_assay_tissue_heatmap.png")

# --- Plot 7: Cell type name similarity heatmap ---
fig, ax = plt.subplots(figsize=(14, 10))
sim_data = []
for et in entex_terms:
    row = []
    for bo in bio_organs:
        score = SequenceMatcher(None, et.lower(), bo.lower()).ratio()
        row.append(score)
    sim_data.append(row)
sim_df = pd.DataFrame(sim_data, index=entex_terms, columns=bio_organs)
# Only show organs with max similarity > 0.4
good_cols = sim_df.columns[sim_df.max() > 0.4]
if len(good_cols) > 0:
    sns.heatmap(sim_df[good_cols], cmap='RdYlGn', vmin=0, vmax=1, ax=ax,
                annot=True, fmt='.2f', linewidths=0.5)
    ax.set_title("Name Similarity: ENTEx Biosample Term Name vs Biosamples.csv Organ/Tissue")
else:
    ax.text(0.5, 0.5, 'No columns with similarity > 0.4', ha='center', va='center')
savefig("07_celltype_similarity_heatmap.png")

# --- Plot 8: Donor × Tissue matrix ---
fig, ax = plt.subplots(figsize=(10, 6))
donor_tissue = meta.groupby(['Donor(s)', 'Biosample term name']).size().reset_index(name='count')
dt_pivot = donor_tissue.pivot(index='Donor(s)', columns='Biosample term name', values='count').fillna(0).astype(int)
# Shorten donor labels
dt_pivot.index = [d.split('/')[-2] if '/' in d else d for d in dt_pivot.index]
sns.heatmap(dt_pivot, annot=True, fmt='d', cmap='Blues', ax=ax, linewidths=0.5)
ax.set_title("ENTEx: File Count by Donor × Biosample Term Name")
savefig("08_donor_tissue_heatmap.png")

# --- Plot 9: Biosamples.csv organ distribution ---
fig, ax = plt.subplots(figsize=(10, 8))
organ_counts = bio['Organ/Tissue'].value_counts().head(30)
organ_counts.plot(kind='barh', ax=ax, color=sns.color_palette("Spectral", len(organ_counts)))
ax.set_title("Top 30 Organ/Tissue Categories in Biosamples.csv")
ax.set_xlabel("Number of Biosamples")
savefig("09_biosamples_organ_distribution.png")

# --- Plot 10: Mapping coverage summary ---
fig, ax = plt.subplots(figsize=(8, 5))
mapping_names = [
    'ENCFF\n(files↔meta)',
    'ENCSR\n(meta↔bio)',
    'ENCFF\n(meta↔bio)',
    'Cell type\n(exact)',
    'Cell type\n(fuzzy→tissue)',
    'Cell type\n(fuzzy→organ)'
]
mapping_vals = [
    len(overlap_file),
    len(overlap_exp),
    len(overlap_meta_bio),
    len(exact_tissue),
    len(fuzzy_matches),
    len(organ_matches)
]
colors = sns.color_palette("viridis", len(mapping_names))
bars = ax.bar(mapping_names, mapping_vals, color=colors)
ax.set_title("Mapping Coverage Summary Across All Linkage Types")
ax.set_ylabel("Number of Matches")
for bar, v in zip(bars, mapping_vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + max(mapping_vals)*0.01, str(v),
            ha='center', fontsize=9)
savefig("10_mapping_coverage_summary.png")

# --- Plot 11: Per-assay experiment overlap between ENTEx and Biosamples ---
fig, ax = plt.subplots(figsize=(8, 5))
assay_map = {
    'DNase-seq': 'DNase Exp. ID',
    'ATAC-seq': 'ATAC Exp. ID',
    'H3K4me3': 'H3K4me3 Exp. ID',
    'H3K27ac': 'H3K27ac Exp. ID',
    'CTCF': 'CTCF Exp. ID'
}
# Map ENTEx assay names to potential bio column matches
entex_assay_exp = {}
for assay in meta['Assay'].unique():
    entex_assay_exp[assay] = set(meta[meta['Assay'] == assay]['Experiment accession'].dropna().unique())

per_assay_overlap = {}
for entex_assay, bio_col in [
    ('DNase-seq', 'DNase Exp. ID'),
    ('ATAC-seq', 'ATAC Exp. ID'),
    ('ChIP-seq', 'H3K4me3 Exp. ID'),  # ChIP-seq covers multiple histone marks
    ('ChIP-seq', 'H3K27ac Exp. ID'),
    ('ChIP-seq', 'CTCF Exp. ID'),
]:
    if entex_assay in entex_assay_exp and bio_col in bio.columns:
        bio_ids = set(bio[bio_col].dropna().unique())
        entex_ids = entex_assay_exp.get(entex_assay, set())
        overlap = entex_ids & bio_ids
        label = f"{entex_assay}↔{bio_col.replace(' Exp. ID','')}"
        per_assay_overlap[label] = len(overlap)

if per_assay_overlap:
    ax.bar(per_assay_overlap.keys(), per_assay_overlap.values(), color=sns.color_palette("Set2"))
    ax.set_title("Per-Assay Experiment ID Overlap (ENTEx ↔ Biosamples)")
    ax.set_ylabel("Shared Experiment IDs")
    plt.xticks(rotation=30, ha='right')
    for i, (k, v) in enumerate(per_assay_overlap.items()):
        ax.text(i, v + 0.5, str(v), ha='center', fontsize=9)
savefig("11_per_assay_experiment_overlap.png")

# ============================================================
# FINAL SUMMARY
# ============================================================
report("\n" + "=" * 80)
report("FINAL SUMMARY & RECOMMENDATIONS")
report("=" * 80)

report(f"""
Key Findings:
1. entex-files.txt contains {len(entex_file_ids_set)} unique file IDs; {len(overlap_file)} ({len(overlap_file)/max(len(entex_file_ids_set),1)*100:.1f}%)
   map to entex-metadata.tsv via File accession.

2. entex-metadata.tsv ↔ Biosamples.csv can be linked via:
   a) Experiment IDs (ENCSR*): {len(overlap_exp)} shared IDs
   b) File IDs (ENCFF*): {len(overlap_meta_bio)} shared IDs
   c) Cell type names: {len(exact_tissue)} exact + {len(fuzzy_matches)} fuzzy matches

3. The strongest mapping key is Experiment accession / Exp. ID (ENCSR* identifiers),
   providing the most reliable cross-reference between tables.

4. Cell type name matching requires normalization because:
   - ENTEx uses short names (e.g., 'ovary', 'thyroid gland')
   - Biosamples.csv uses descriptive names (e.g., 'Ovary, female adult (53 years)')
   - Biosamples.csv Organ/Tissue provides broader category matches

5. The Biosample term id (UBERON:*/CL:* ontology) in entex-metadata provides
   standardized identifiers but is NOT present in Biosamples.csv.

6. All three tables share ENCODE file/experiment accession IDs as the primary
   linkage mechanism. URL patterns in both entex-files.txt and Biosamples.csv
   contain these IDs embedded in download links.
""")

save_report()
print(f"\nAll outputs saved to: {OUT_DIR}/")
