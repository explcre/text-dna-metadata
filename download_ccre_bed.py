#!/usr/bin/env python3
"""Download all bed-format cCRE files listed in entex-metadata.tsv."""

import csv
import os
import subprocess
import sys
from collections import Counter

METADATA = "entex-metadata.tsv"
OUTDIR = "ccre_bed_files"


def parse_bed_entries(metadata_path):
    """Parse metadata TSV and return rows where File type == 'bed'."""
    entries = []
    with open(metadata_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["File type"] == "bed":
                entries.append(row)
    return entries


def download_file(url, dest):
    """Download a file using wget. Returns True on success."""
    try:
        subprocess.run(
            ["wget", "-q", "-O", dest, url],
            check=True, timeout=300,
        )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def main():
    if not os.path.exists(METADATA):
        sys.exit(f"Error: {METADATA} not found in current directory.")

    entries = parse_bed_entries(METADATA)
    total = len(entries)
    print(f"Found {total} bed files in {METADATA}")

    # Print metadata statistics
    assays = Counter(e["Assay"] for e in entries)
    biosamples = Counter(e["Biosample term name"] for e in entries)
    output_types = Counter(e["Output type"] for e in entries)
    assemblies = Counter(e["File assembly"] for e in entries)

    print(f"\n--- Metadata Statistics ---")
    print(f"Total bed files: {total}")
    print(f"Unique assays: {len(assays)}")
    for a, c in assays.most_common():
        print(f"  {a}: {c}")
    print(f"Unique biosamples: {len(biosamples)}")
    for b, c in biosamples.most_common(10):
        print(f"  {b}: {c}")
    if len(biosamples) > 10:
        print(f"  ... and {len(biosamples) - 10} more")
    print(f"Output types:")
    for o, c in output_types.most_common():
        print(f"  {o}: {c}")
    print(f"Assemblies: {dict(assemblies)}")

    # Download
    os.makedirs(OUTDIR, exist_ok=True)
    success = 0
    failed = []
    total_size = 0

    for i, entry in enumerate(entries, 1):
        accession = entry["File accession"]
        url = entry["File download URL"]
        filename = url.split("/")[-1]
        dest = os.path.join(OUTDIR, filename)

        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"[{i}/{total}] {filename} already exists, skipping")
            success += 1
            total_size += os.path.getsize(dest)
            continue

        print(f"[{i}/{total}] Downloading {filename} ...", end=" ", flush=True)
        if download_file(url, dest):
            fsize = os.path.getsize(dest)
            total_size += fsize
            success += 1
            print(f"OK ({fsize:,} bytes)")
        else:
            failed.append(accession)
            print("FAILED")

    # Summary
    print(f"\n--- Download Summary ---")
    print(f"Successful: {success}/{total}")
    print(f"Failed: {len(failed)}")
    if failed:
        print(f"Failed accessions: {', '.join(failed)}")
    print(f"Total downloaded size: {total_size / 1e9:.2f} GB")
    print(f"Files saved to: {OUTDIR}/")


if __name__ == "__main__":
    main()
