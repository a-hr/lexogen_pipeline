#!/usr/bin/env python

import argparse
import os
import re
from pathlib import Path

import pandas as pd


def parse() -> tuple:
    parser = argparse.ArgumentParser(
        description="Demultiplex .fastq files given the sample and experiment barcodes."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Experiment and group barcode .csv files' directory",
        required=True,
    )
    parser.add_argument(
        "-p", "--path", type=str, help=".fastq files directory", required=True
    )
    parser.add_argument(
        "-s",
        "--suffixes",
        help="Suffixes of R1 and R2 strands",
        nargs=2,
        default=["_R1", "_R2"],
    )

    _args = parser.parse_args()
    return _args.input, _args.path, _args.suffixes


def demultiplex(_input: Path, fastq: Path, _suffixes: tuple = ("_R1", "_R2")) -> None:
    fp_table = pd.read_csv(
        Path(_input) / "fiveprime.csv",
        header=None,
        names=["fastq", "exp"],
        sep=";",
    ).dropna(axis=0, inplace=False)

    tp_table = pd.read_csv(
        Path(_input) / "threeprime.csv",
        header=None,
        names=["id", "barcode", "sub_exp"],
        sep=";",
    ).dropna(axis=0, inplace=False)

    bcs: list = [
        f"{(row[1]['sub_exp']).strip()}=^{(row[1]['barcode']).strip()}"
        for row in tp_table[["sub_exp", "barcode"]].iterrows()
    ]
    bcs: str = " -g ".join(bcs)  # contains the barcodes and their id, concatenated with -g

    fq1 = [str(fq) for fq in Path(fastq).iterdir() if _suffixes[0] in str(fq)]
    fq2 = [str(fq) for fq in Path(fastq).iterdir() if _suffixes[1] in str(fq)]

    # sort both lists so that files are paired (R1/R2)
    fq1.sort()
    fq2.sort()

    if not (fq1 and fq2):
        print("No .fastq files provided or wrong suffixes, check the path")
        quit()

    exp_files_list = _get_exps(fp_table, fq1, fq2)

    for i, data in enumerate(exp_files_list.items()):
        exp, r = data
        cmd = f"cutadapt -e 0.2 --no-indels --cores 0 -g {bcs} -p {exp}_{{name}}_R1.fastq.gz -o {exp}_{{name}}_R2" \
              f".fastq.gz {r[1]} {r[0]}> demultiplex.{i}.log"
        os.system(cmd)


def _get_exps(fp_table, _fq1, _fq2) -> dict:
    """Given a list of R1 or R2 files, gets their corresponding experiment id from the table.
    The response has the following format:
        {
            exp: (R1, R2),
            ...
        }
    """
    fq1, fq2 = _fq1.copy(), _fq2.copy()
    response = dict()

    for _, row in fp_table.iterrows():
        for index, file in enumerate(fq1):
            r = re.search(f"({row.fastq})", str(file))
            if r:  # match found
                response[row.exp] = (fq1.pop(index), fq2.pop(index))
                break
    return response


if __name__ == "__main__":
    args = parse()
    demultiplex(*args)  # type: ignore
