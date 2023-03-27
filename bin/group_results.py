#!/usr/bin/env python
from pathlib import Path

import click
import pandas as pd


def _create_weights(samples: list) -> dict:
    d = {}
    w = 1
    for s in samples:
        d[s] = w
        w += 1
    return d


def _classify(sample: str, ref: dict) -> str:
    ref = ref.keys()
    sample = sample.upper()

    for r in ref:
        if r in sample:
            return r


def _sort_columns(col_name: str, mouses: list, tissues: list):

    if not "_" in col_name:
        # ctrl sample
        total_weight = 999

    else:
        # Assign a weight to each mouse type and tissue
        mouse_weights = _create_weights(mouses)
        tissue_weights = _create_weights(tissues)

        mouse_type = _classify(col_name, mouse_weights)
        tissue = _classify(col_name, tissue_weights)

        # Calculate the total weight of each column name
        total_weight = mouse_weights[mouse_type] * 10 + tissue_weights[tissue]

    return total_weight


def group_csv(path: str, mouses: list, tissues: list) -> pd.DataFrame:
    """Read a csv and reorder the columns in alphabetical order."""
    df = pd.read_csv(path, index_col=0, sep="\t")
    info = df[["Chr", "Start", "End", "Length", "Strand"]]
    data = df.drop(["Chr", "Start", "End", "Length", "Strand"], axis=1)

    data = data.reindex(
        sorted(data.columns, key=lambda col: _sort_columns(col, mouses, tissues)),
        axis=1,
    )

    df = pd.concat([info, data], axis=1)
    return df


@click.command()
@click.option("--output", "-o", default="results", help="Path to the output dir.")
@click.option("--sep", "-s", default="\t", help="Separator used in the output file.")
@click.option("--mouse", "-m", multiple=True, help="Mouse types to group.")
@click.option("--tissue", "-t", multiple=True, help="Tissue types to group.")
@click.argument("paths", nargs=-1)
def main(output, sep, mouse, tissue, **paths):
    paths = paths["paths"]
    for path in paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"File {path} was not found")
    else:
        print(f"Grouping the following files:\n{paths}")

    outdir = Path.cwd()
    if not (rdir := (Path(outdir) / output)).exists():
        rdir.mkdir(parents=True)
        print(f"Output dir {output} not found.")
        print("Creating the results directory...")

    for path in paths:
        outname = Path(path).name
        df = group_csv(path, mouse, tissue)
        df.to_csv(f"{rdir}/{outname}", sep=sep)
    print("Done.")


if __name__ == "__main__":
    """example command:
    python group_results.py files/*.csv --params params.yaml
    """
    main()
