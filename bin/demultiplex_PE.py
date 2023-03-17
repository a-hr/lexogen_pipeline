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
        sys.stderr.write("InputError: Barcode csv does not exist\n")
        quit()

    if not (fs_path := Path(fs)).exists():
        sys.stderr.write("InputError: .fastq files directory does not exist\n")
        quit()

    bc_df = pd.read_csv(bc_path, sep=";").dropna(axis=0, inplace=False)
    bc_df.Sample.replace(" ", "_", inplace=True, regex=True)

    try:
        f1 = next(fs_path.glob(f"*{s[0]}.fastq.gz"))
    except StopIteration:
        sys.stderr.write("InputError: No R1 .fastq file found\n")
        quit()
    
    if not (f2 := Path(str(f1).replace(s[0], s[1]))).exists():
        sys.stderr.write("InputError: No R2 .fastq file found\n")
        quit()

    return bc_df, f1, f2

def demultiplex(bc_csv: Path, f1: Path, f2: Path) -> None:
    bcs: list = [
        f"{(row[1]['Sample']).strip()}=^{(row[1]['Sequence']).strip()}"
        for row in bc_csv.iterrows()
    ]
    bcs = " -g ".join(bcs)  # contains the barcodes and their id, concatenated with -g

    cmd = f"""cutadapt -e 0.2 --no-indels -j 0 -g {bcs} -p {{name}}_R1.fastq.gz -o {{name}}_R2.fastq.gz {f2} {f1}"""
    os.system(cmd)


if __name__ == "__main__":
    args = parse()
    bc_df, f1, f2 = validate_inputs(*args)
    demultiplex(bc_df, f1, f2)
