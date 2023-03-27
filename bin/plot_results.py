#!/usr/bin/env python
import re
from pathlib import Path

import click
import pandas as pd
from plotly.express import box


# -------------- HELPER FUNCTIONS -------------- #
def _classify(sample: str, ref: list) -> str:
    sample = sample.upper()
    for r in ref:
        if r in sample:
            return r


def _parse_ctrl(df, label, match, control, groups) -> pd.DataFrame:
    if groups == 3:
        i7, sample = match.groups()
    elif groups == 2:
        sample = label

    if groups == 3:
        df.loc[label, "lib"] = i7
    df.loc[label, "tissue"] = "control"
    df.loc[label, "sample"] = _classify(sample, control)

    return df


def _parse(df, label, samples, control, groups) -> pd.DataFrame:
    # ----- PARSING -----
    if groups == 3:
        # parse data if 3 groups expected
        if match := re.match(r"^(.*)_(.*)_(.*)$", label):
            i7, sample, tissue = match.groups()
        else:
            # if ctrl sample, let aux func handle parsing
            match = re.match(r"^(.*)_(.*)$", label)
            return _parse_ctrl(df, label, match, control, groups)

    elif groups == 2:
        # parse data if 2 groups expected
        if match := re.match(r"^(.*)_(.*)$", label):
            sample, tissue = match.groups()
        else:
            # if ctrl sample, let aux func handle parsing
            match = re.match(r"^(.*)_(.*)$", label)
            return _parse_ctrl(df, label, match, control, groups)

    # ----- ASSIGNING -----
    if not (s := _classify(sample, samples)):
        return df

    if groups == 3:
        df.loc[label, "lib"] = i7
    df.loc[label, "tissue"] = tissue.upper()
    df.loc[label, "sample"] = s

    return df


def preprocess_data(df, samples, target, control, groups):
    # prepare the df to contain the required data and drop useless cols
    drop_cols = ["Chr", "Start", "End", "Strand", "Length"]

    df = df.drop(drop_cols, axis=1, inplace=False)
    df = df.set_index("Geneid").T

    if groups == 3:
        df["lib"] = ""
    df["sample"] = ""
    df["tissue"] = ""
    df["target"] = ""

    # parse the index (filename) to extract the required data
    for label, _ in df.iterrows():
        df = _parse(df, label, samples, control, groups)

    # for each row in the df, separate the row in 3 df depending on its target
    final_df = pd.DataFrame(columns=["lib", "sample", "tissue", "target", "counts"])
    for t in target:
        rest = [r for r in target if r != t]
        # create a copy of df without the columns in rest list
        df_t = df.drop(rest, axis=1, inplace=False)
        df_t["target"] = t
        df_t.rename(columns={t: "counts"}, inplace=True)
        # append the new df to the final df
        final_df = pd.concat([final_df, df_t])
    return final_df


def plot_data(df, title):
    # dont show control data
    # df = df[df["tissue"] != "control"]
    # add 1 to counts to avoid log(0)
    df = df.assign(counts=df["counts"] + 1)
    # add facet row to show the data for each target
    fig = box(
        df,
        x="sample",
        y="counts",
        log_y=True,
        color="tissue",
        points="all",
        facet_col="target",
    )
    # increase the point size
    fig.update_traces(marker_size=3)

    # increase legend font size
    fig.update_layout(legend=dict(font=dict(size=18)))

    # add title
    fig.update_layout(title_text=title)

    # remove the target= part from the facet row title
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    # change axis label
    fig.update_yaxes(title_text="log(counts + 1)", row=1, col=1)

    return fig


def export_plot(fig, outdir, title):
    if not (p := Path(outdir)).resolve().exists():
        print("Output dir '{outdir}' not found. Creating...")
        p.mkdir(parents=True)

    fig.write_image(f"{p}/{title}.png", width=1920, height=1080, scale=6)
    fig.write_html(f"{p}/{title}.html")


# -------------- MAIN -------------- #
@click.command()
@click.option("--lib", "-l", help="Lib name to be added to the libs to group by.")
@click.option("--sample", "-s", multiple=True, help="Sample name.")
@click.option("--target", "-t", multiple=True, help="Must appear as in the .tsv file.")
@click.option("--control", "-c", multiple=True, help="Control sample name.")
@click.option("--pattern", "-p", default=3, help="Pattern groups")
@click.option("--outdir", "-o", default="results", help="Path to the output dir.")
@click.argument("files", nargs=-1, type=click.Path(exists=True))
def main(lib, sample, target, control, outdir, pattern, **files):
    # todo: add an extra plot by library. currently ignored argument (usually one pipeline run by lib).

    # -------- INPUT CHECKING --------
    files = files["files"]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(f"File {file} was not found")
    else:
        print(f"Plotting the following files:\n{files}")

    if isinstance(lib, str):
        lib = [
            lib,
        ]

    if isinstance(control, str):
        control = [
            control,
        ]

    if isinstance(sample, str):
        sample = [
            sample,
        ]

    # --------  PROCESSING --------
    for f in files:
        # read the data
        counts = pd.read_csv(f, sep="\t")
        counts_processed = preprocess_data(counts, sample, target, control, pattern)

        # plot the data
        title = Path(f).stem
        p = plot_data(counts_processed, title)

        # export
        export_plot(p, outdir, title)


if __name__ == "__main__":
    """Exaple command:
        python plot_results.py inputs/* -l lib1 -t hHTT_PA1_v1 -t hHTT_PA2_v1 -t hHTT_cPA1_v1 -t hHTT_cPA2_v1 -s WT -s YAC -s KI -c NTC -p 2
    """
    main()
