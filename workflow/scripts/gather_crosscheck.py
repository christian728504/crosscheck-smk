"""Gather step: concatenate a job's per-chunk CrosscheckFingerprints TSVs into one
crosscheck_{job}.tsv.

Cross jobs (LEFT != RIGHT modality): each LEFT sample is in exactly one chunk, so
the LEFT x RIGHT pairs partition cleanly -- a plain vertical concat.

Self jobs (LEFT == RIGHT): the chunk's samples appear in both --INPUT and
--SECOND_INPUT, so Picard emits self-pairs (a x a) and both orientations of each
unordered pair. When snakemake.params.dedup is true we drop self-pairs and keep one
row per unordered pair.

Picard metrics TSVs carry '#'/'##' comment lines AND a blank line before the table,
which trips polars' schema inference; strip both before parsing (same idiom as
src/detect_sample_swap/plot.py). The output is a clean header+rows TSV that
aggregate.py reads downstream. Run via Snakemake's ``script:`` directive.
"""

from pathlib import Path

import polars as pl


def read_crosscheck(path: str) -> pl.DataFrame:
    raw = "".join(
        line
        for line in Path(path).open()
        if line.strip() and not line.startswith("#")
    ).encode("utf-8")
    return pl.read_csv(raw, separator="\t", infer_schema_length=None)


df = pl.concat(
    [read_crosscheck(p) for p in snakemake.input.chunks],  # noqa: F821 - injected
    how="vertical",
)

if snakemake.params.dedup:  # noqa: F821 - self comparison
    df = (
        df.filter(pl.col("LEFT_GROUP_VALUE") != pl.col("RIGHT_GROUP_VALUE"))
        .with_columns(
            pl.min_horizontal("LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE").alias("_a"),
            pl.max_horizontal("LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE").alias("_b"),
        )
        .unique(subset=["_a", "_b"], keep="first")
        .drop("_a", "_b")
    )

df.write_csv(snakemake.output.tsv, separator="\t")  # noqa: F821
