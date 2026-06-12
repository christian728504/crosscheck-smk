"""Pre-flight: build the dataset-sample-participant manifest for crosscheck-smk.

Standalone (run BEFORE snakemake) because the sample set must exist when snakemake
builds the DAG. Reproduces investigate-father-son/construct-manifest.py:

  1. Read the GCP-Upload dataset.*.tsv and sample.*.tsv tables.
  2. Strip whitespace from column names and string cells; dedup.
  3. Drop omics modalities (lipidomics/metabolomics/proteomics/exposomics).
  4. Left-join dataset.entity_id == sample."entity:sample_id"; drop rows missing
     "entity:dataset_id" or participant_id; sort by "entity:dataset_id".
  5. Resolve each row's alignment path by modality and keep only existing files.
  6. Write a parquet with columns:
     entity:dataset_id, data_modality, participant_id, file, file_exists.

Needs polars + pyyaml (both in snakemake-base), e.g.
`conda run -n snakemake-base python workflow/scripts/construct_manifest.py`.
"""

import argparse
from pathlib import Path

import polars as pl
import yaml

# Modality -> alignment path template. Each {} is filled with entity:dataset_id
# (passed twice to pl.format). Mirrors the notebook; WGS is CRAM, the rest BAM.
PATH_TEMPLATES = {
    "WGS": "/data/projects/mohd/data/Molecular/0_WGS/{}/{}_align_GRCh38_v0.cram",
    "WGBS": "/data/projects/mohd/data/Molecular/1_WGBS/{}/{}_align_GRCh38_v0.bam",
    "ATAC-seq": "/data/projects/mohd/data/Molecular/2_ATAC/{}/{}_unfiltered_GRCh38_v0.bam",
    "RNA-seq": "/data/projects/mohd/data/Molecular/3_RNA/{}/{}_align_GRCh38_v0.bam",
}

DROP_MODALITIES = ["lipidomics", "metabolomics", "proteomics", "Exposomics", "exposomics"]


def _clean(df: pl.DataFrame) -> pl.DataFrame:
    """Strip column names and string cells, then drop duplicate rows."""
    renamed = {c: c.strip() for c in df.columns}
    return (
        df.rename(renamed)
        .with_columns(pl.selectors.string().str.strip_chars(" "))
        .unique()
    )


def build_manifest(data_model_dir: Path, exclude=()) -> pl.DataFrame:
    dataset_df = _clean(
        pl.read_csv(str(data_model_dir / "dataset.*.tsv"), separator="\t")
    ).filter(pl.col("data_modality").is_in(DROP_MODALITIES).not_())

    sample_df = _clean(
        pl.read_csv(str(data_model_dir / "sample.*.tsv"), separator="\t")
    )

    join_df = (
        dataset_df.join(
            sample_df,
            how="left",
            left_on="entity_id",
            right_on="entity:sample_id",
        )
        .drop_nulls(subset=["entity:dataset_id", "participant_id"])
        .sort("entity:dataset_id")
    )

    # Build the alignment path per row from its modality template.
    file_expr = pl.when(pl.col("data_modality") == "WGS").then(
        pl.format(PATH_TEMPLATES["WGS"], "entity:dataset_id", "entity:dataset_id")
    )
    for modality in ("WGBS", "ATAC-seq"):
        file_expr = file_expr.when(pl.col("data_modality") == modality).then(
            pl.format(PATH_TEMPLATES[modality], "entity:dataset_id", "entity:dataset_id")
        )
    file_expr = file_expr.otherwise(  # RNA-seq
        pl.format(PATH_TEMPLATES["RNA-seq"], "entity:dataset_id", "entity:dataset_id")
    )

    file_df = join_df.with_columns(file_expr.alias("file"))
    file_df = file_df.with_columns(
        pl.col("file")
        .map_elements(lambda p: Path(p).exists(), return_dtype=pl.Boolean)
        .alias("file_exists")
    ).filter("file_exists")

    if exclude:
        file_df = file_df.filter(pl.col("entity:dataset_id").is_in(list(exclude)).not_())

    return file_df.select(
        "entity:dataset_id", "data_modality", "participant_id", "file", "file_exists"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Pipeline config (for data_model_dir and manifest defaults).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output parquet path (default: config['manifest']).",
    )
    parser.add_argument(
        "--data-model-dir",
        default=None,
        help="Override config['data_model_dir'].",
    )
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    data_model_dir = Path(args.data_model_dir or config["data_model_dir"])
    out = Path(args.out or config["manifest"])
    exclude = config.get("exclude_samples") or []

    manifest = build_manifest(data_model_dir, exclude)
    out.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_parquet(out)

    if exclude:
        print(f"excluded {len(exclude)} sample(s): {exclude}")
    print(f"wrote {manifest.height} samples -> {out}")
    print(manifest["data_modality"].value_counts())


if __name__ == "__main__":
    main()
