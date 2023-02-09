#!/usr/bin/env python

import os
import sys
import argparse

from pathlib import Path

def parse() -> dict:
    parser = argparse.ArgumentParser(
        description="STAR aligner wrapper."
    )
    parser.add_argument(
        "-t", "--task", type=str, help="run mode (index/align)", choices=["index", "align"], required=True
    )
    # * index generation args
    parser.add_argument("-I", "--index", type=str, help="index folder") # req. for align too
    parser.add_argument("-g", "--genome", type=str, help="ref. genome directory")
    parser.add_argument("-L", "--len", type=int, help="sjdbOverhang length")

    # * alignment args
    parser.add_argument("-f", "--file", type=str, help="file to align")
    parser.add_argument("-o", "--output", type=str, help="output dir", default="")
    parser.add_argument("-s", "--setts", type=str, help="file containing extra args in .json format", default="")

    args = parser.parse_args()
    if args.task == "index":
        return dict(
            task=args.task, 
            index=args.index,
            genome=args.genome,
            length=args.len
            ) 
    return dict(
            task=args.task, 
            index=args.index,
            file=args.file,
            out=args.output,
            setts=args.setts
            )

def runSTAR(**kwargs) -> None:
    if kwargs["task"] == "index":
        cmd = _STAR_index(**kwargs)
        sys.stderr.write("index task")
    else:
        cmd = _STAR_align(**kwargs)
        sys.stderr.write("align task")
    sys.stderr.write(cmd)
    os.system(cmd)

def _STAR_index(**kwargs) -> str:
    gDir = kwargs.get("index")
    fa_file = Path(kwargs.get("genome")).glob("*.fa")[0] # type: ignore
    annot = Path(kwargs.get("genome")).glob("*.gtf")[0] # type: ignore
    length = kwargs.get("length")

    return f"""STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir {gDir} \
    --genomeFastaFiles {str(fa_file)} \
    --sjdbGTFfile {str(annot)} \
    --sjdbOverhang {length}
    """

def _STAR_align(**kwargs) -> str:
    """Performs STAR alignment with shared memory to allow parallel runs. STAR will check if genome
     is loaded, and load it only if not. Every time a job finishes, it will check if any other job 
     is running. If not the case, it will remove the index from RAM. If an error crashes the program,
     shared memory is not freed and the following command must be run:
        STAR --genomeLoad Remove --genomeDir $STARINDEX
    """
    gDir = kwargs.get("index")
    file = kwargs.get("file")

    out_file: str = file # type: ignore
    for s in (".gz", ".fq", ".fastq"):
        out_file = out_file.replace(s, "")

    if (d:=kwargs.get("out")) != "":
        out_file = str(Path(d) / out_file) # type: ignore
    #     --genomeLoad LoadAndRemove \
    return f"""STAR --runThreadN 8 \
    --genomeDir {gDir} \
    --readFilesIn {file} \
    --readFilesCommand zcat \
    --outFileNamePrefix {out_file} \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --limitBAMsortRAM 10000000000 \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.04 \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 \
    --alignEndsType Local \
    --outSAMattributes Standard 
    """

def main() -> None:
    kwargs = parse()
    runSTAR(**kwargs)

if __name__ == "__main__":
    main()