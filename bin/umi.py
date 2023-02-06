#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path

def parse() -> tuple:
    parser = argparse.ArgumentParser(
        description="Handle UMI sequences on .fastq files."
    )
    parser.add_argument(
        "-m", "--mode", type=str, help="run mode (extract/dedup)", choices=["extract", "dedup"], required=True
    )
    parser.add_argument(
        "-I", "--input", type=str, help="input .fastq file", required=True
    )
    parser.add_argument(
        "-P", "--pattern", type=str, default="--bc-pattern NNNNNN", help="barcode pattern to extract"
    )

    args = parser.parse_args()
    return (args.mode, args.input, args.pattern)

def umi_run(mode: str, file: str, pattern: str) -> None:
    p = pattern if mode == "extract" else ""
    outfile = file
    exts = Path(outfile).suffixes
    for ext in exts:
        outfile = outfile.replace(ext, "")
    exts = "".join(exts)
    outfile = f"{outfile}_UMI{exts}"

    cmd = f"umi_tools {mode} {p} -I {file} -S {outfile}"
    sys.stderr.write(cmd)
    os.system(cmd)

if __name__ == "__main__":
    args = parse()
    umi_run(*args) # type: ignore