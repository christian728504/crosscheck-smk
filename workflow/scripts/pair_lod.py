"""Per-pair haplotype-block LOD decomposition for one flagged sample pair.

Calls snp_prioritization.pair_lod on the pair's two fingerprint VCFs and the job's
haplotype map (extraction map == crosscheck map, so block representatives align), and
writes the per-block contributions to a parquet. Σ delta reproduces Picard's
LOD_SCORE for the pair; ranking by delta surfaces the blocks driving a swap call.

Pair identity (left, right, comparison) is encoded in the output path, not the table.
Run via Snakemake's ``script:`` directive (``snakemake`` object injected).
"""

import polars as pl
from snp_prioritization import pair_lod

SCHEMA = {
    "chrom": pl.String,
    "pos": pl.Int64,
    "name": pl.String,
    "maf": pl.Float64,
    "delta": pl.Float64,
}

result = pair_lod(
    snakemake.input.vcf1,  # noqa: F821 - injected by Snakemake
    snakemake.input.vcf2,  # noqa: F821
    snakemake.input.hmap,  # noqa: F821
)

df = pl.DataFrame(
    {
        "chrom": result.chrom.tolist(),
        "pos": result.pos.tolist(),
        "name": result.name.tolist(),
        "maf": result.maf.tolist(),
        "delta": result.delta.tolist(),
    },
    schema=SCHEMA,
)
df.write_parquet(snakemake.output.parquet)  # noqa: F821
