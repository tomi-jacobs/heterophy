#!/usr/bin/env python3
"""
validate_samplesheet.py
────────────────────────
HeteroPhy — Step 0: Validate and normalise the input samplesheet.
Checks that all required columns are present and that all files exist.
"""

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {"sample", "reads_r1", "reads_r2", "transcriptome"}


def parse_args():
    parser = argparse.ArgumentParser(description="Validate HeteroPhy samplesheet")
    parser.add_argument("--input",  required=True, help="Input CSV samplesheet")
    parser.add_argument("--output", required=True, help="Output validated CSV")
    return parser.parse_args()


def main():
    args = parse_args()
    errors = []
    valid_rows = []

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])

        missing_cols = REQUIRED_COLUMNS - columns
        if missing_cols:
            print(f"[HeteroPhy] ERROR: Samplesheet missing required columns: "
                  f"{', '.join(sorted(missing_cols))}")
            print(f"[HeteroPhy] Required columns: {', '.join(sorted(REQUIRED_COLUMNS))}")
            sys.exit(1)

        for i, row in enumerate(reader, start=2):
            row_errors = []

            # Check sample name
            sample = row.get("sample", "").strip()
            if not sample:
                row_errors.append("'sample' is empty")

            # Check files exist
            for col in ["reads_r1", "reads_r2", "transcriptome"]:
                path = row.get(col, "").strip()
                if not path:
                    row_errors.append(f"'{col}' is empty")
                elif not Path(path).exists():
                    row_errors.append(f"File not found: {path}")

            if row_errors:
                errors.append(f"  Row {i} ({sample}): {'; '.join(row_errors)}")
            else:
                valid_rows.append({
                    "sample": sample,
                    "reads_r1": row["reads_r1"].strip(),
                    "reads_r2": row["reads_r2"].strip(),
                    "transcriptome": row["transcriptome"].strip(),
                })

    if errors:
        print(f"[HeteroPhy] Samplesheet validation FAILED with {len(errors)} error(s):")
        for err in errors:
            print(err)
        sys.exit(1)

    # Write validated samplesheet
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "reads_r1", "reads_r2", "transcriptome"])
        writer.writeheader()
        writer.writerows(valid_rows)

    print(f"[HeteroPhy] Samplesheet validated: {len(valid_rows)} sample(s) passed")
    print(f"[HeteroPhy] Validated samplesheet written to: {args.output}")


if __name__ == "__main__":
    main()
