#!/usr/bin/env python3

import argparse
import csv
import glob
import sys
import shutil
import os

def main():
    parser = argparse.ArgumentParser(
        description="Select the appropriate runsheet for a given target region."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing runsheet CSV files")
    parser.add_argument("--target", required=True,
                        help="Target region (e.g. 16S, 18S, ITS)")
    parser.add_argument("--specify-runsheet",
                        help="Specify a runsheet by name when multiple match the same target.",
                        default=None)
    args = parser.parse_args()

    target = args.target.lower()
    runsheet_files = glob.glob(f"{args.input_dir}/*.csv")

    if len(runsheet_files) == 0:
        sys.exit(f"Error: No runsheet CSV files found in {args.input_dir}/.")

    if args.specify_runsheet:
        match = f"{args.input_dir}/{args.specify_runsheet}"
        if match in runsheet_files:
            shutil.copy(match, os.path.basename(match))
            print(f"Using specified runsheet: {args.specify_runsheet}")
            sys.exit(0)
        else:
            sys.exit(
                f"Error: Specified runsheet '{args.specify_runsheet}' not found. "
                f"Available: {[os.path.basename(f) for f in runsheet_files]}"
            )

    if len(runsheet_files) == 1:
        shutil.copy(runsheet_files[0], os.path.basename(runsheet_files[0]))
        print(f"Single runsheet found, using: {os.path.basename(runsheet_files[0])}")
        sys.exit(0)

    # Multiple runsheets — filter by target region
    def get_target_region(filepath):
        with open(filepath, newline='') as f:
            reader = csv.DictReader(f)
            regions = set()
            for row in reader:
                regions.add(row.get("Parameter Value[Library Selection]", "").lower())
            return regions

    matching = []
    for f in runsheet_files:
        try:
            regions = get_target_region(f)
            if len(regions) == 1 and target in regions:
                matching.append(f)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)

    if len(matching) == 1:
        shutil.copy(matching[0], os.path.basename(matching[0]))
        print(f"Using runsheet: {os.path.basename(matching[0])}")
    elif len(matching) > 1:
        sys.exit(
            f"Multiple runsheets match target region '{target}': "
            f"{[os.path.basename(f) for f in matching]}.\n"
            f"Use --specify-runsheet to select one."
        )
    else:
        sys.exit(
            f"No runsheet matches target region '{target}'. "
            f"Available: {[os.path.basename(f) for f in runsheet_files]}"
        )

if __name__ == "__main__":
    main()