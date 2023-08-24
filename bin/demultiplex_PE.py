#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def parse() -> tuple:
    parser = argparse.ArgumentParser(
        description="Demultiplex .fastq files given the sample and experiment barcodes."
    )
    parser.add_argument("-b", type=str, help="Barcode table (csv) directory path")
    parser.add_argument("-f", type=str, help=".fastq files directory")
    parser.add_argument("-s", nargs=2, default=["_R1", "_R2"])

    _args = parser.parse_args()
    return _args.b, _args.f, _args.s

def validate_inputs(bc: str, fs: str, s: list) -> tuple:
    if not (bc_path := next(Path(bc).glob("threeprime.csv"))).exists():
        sys.stderr.write("InputError: Barcode (3') csv does not exist\n")
        quit()

    if not (exp_path := next(Path(bc).glob("fiveprime.csv"))).exists():
        sys.stderr.write("InputError: Experiment (5') csv does not exist\n")
        quit()

    if not (fs_path := Path(fs)).exists():
        sys.stderr.write("InputError: .fastq files directory does not exist\n")
        quit()

    bc_df = pd.read_csv(bc_path, sep=";").dropna(axis=0, inplace=False)
    bc_df.Sample.replace(" ", "_", inplace=True, regex=True)
    bc_df.Sample.replace("/", "-", inplace=True, regex=True)

    exp_df = pd.read_csv(exp_path, sep=";").dropna(axis=0, inplace=False)

    try:
        f1 = next(fs_path.glob(f"*{s[0]}.fastq.gz"))
    except StopIteration:
        sys.stderr.write("InputError: No R1 .fastq file found\n")
        quit()
    
    # replace the last occurence of the suffix to avoid replacing the suffix in the filename
    if not (f2 := Path(s[1].join(str(f1).rsplit(s[0], 1)))).exists():
        sys.stderr.write("InputError: No R2 .fastq file found\n")
        quit()

    return bc_df, exp_df, f1, f2, s

def demultiplex(bc_csv: Path, exp_csv: Path, f1: Path, f2: Path, suffix: list) -> None:
    bcs: list = [
        f"{(row[1]['Sample']).strip()}=^{(row[1]['Sequence']).strip()}"
        for row in bc_csv.iterrows()
    ]
    bcs = " -g ".join(bcs)  # contains the barcodes and their id, concatenated with -g

    sample_id = f1.stem.split(".")[0].removesuffix(suffix[0]).removeprefix("trim_")
    exp = exp_csv.loc[exp_csv["File"] == sample_id].iloc[0, 1]
    
    cmd = f"""cutadapt -e 0 --no-indels -j 0 -m 25 -g {bcs} --json={exp}_bc.cutadapt.json -p {exp}_{{name}}_R1.fastq.gz -o {exp}_{{name}}_R2.fastq.gz {f2} {f1}"""
    os.system(cmd)


if __name__ == "__main__":
    args = parse()
    bc_df, exp_df, f1, f2, suffix = validate_inputs(*args)
    demultiplex(bc_df, exp_df, f1, f2, suffix)
