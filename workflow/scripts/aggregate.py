"""Concatenate per-job CrosscheckFingerprints TSVs and rank cross-participant pairs.

Writes:
  - combined.parquet: the kept columns from every job plus a COMPARISON column.
  - cross_participant.tsv: pairs whose two samples belong to different individuals
    (LEFT_GROUP_VALUE != RIGHT_GROUP_VALUE), sorted by LOD_SCORE descending --
    candidate relatives. Each unordered pair appears twice; that is expected.

Run via Snakemake's ``script:`` directive (``snakemake`` object injected).
"""

from pathlib import Path

import polars as pl

KEEP = [
    "LEFT_GROUP_VALUE",
    "RIGHT_GROUP_VALUE",
    "RESULT",
    "DATA_TYPE",
    "LOD_SCORE",
    "LEFT_FILE",
    "RIGHT_FILE",
]


def comparison_name(path: str) -> str:
    """job name from results/crosscheck/crosscheck_<job>.tsv."""
    return Path(path).stem.removeprefix("crosscheck_")


frames = [
    pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=None)
    .select(KEEP)
    .with_columns(pl.lit(comparison_name(path)).alias("COMPARISON"))
    for path in snakemake.input.tsvs  # noqa: F821 - injected by Snakemake
]

combined = pl.concat(frames, how="vertical")
combined.write_parquet(snakemake.output.combined)  # noqa: F821

(
    combined.filter(pl.col("LEFT_GROUP_VALUE") != pl.col("RIGHT_GROUP_VALUE"))
    .sort("LOD_SCORE", descending=True)
    .write_csv(snakemake.output.cross, separator="\t")  # noqa: F821
)
