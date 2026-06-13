# Per-block LOD decomposition for the flagged sample pairs (checkpoint fan-out).
#
# `cases` is a checkpoint: the set of flagged pairs is only known once cases.tsv
# exists, so the number of per-pair parquets is dynamic. For each UNEXPECTED_* row in
# cases.tsv we run snp_prioritization.pair_lod on that pair's two fingerprint VCFs and
# the job's haplotype map, writing the per-block deltas to
# results/pair_lod/{job}/{left}~{right}.parquet (pair identity is in the path).


wildcard_constraints:
    # dataset_ids never contain a path separator or the '~' pair delimiter.
    left=r"[^/~]+",
    right=r"[^/~]+",


def pair_lod_inputs(wildcards):
    # Both fingerprints were extracted with the job's map (extraction map == crosscheck
    # map), so the block representatives line up with that same map for pair_lod.
    m = job_map(wildcards.job)
    vcf1, vcf2 = fp(wildcards.left, m), fp(wildcards.right, m)
    return {
        "vcf1": vcf1,
        "vcf2": vcf2,
        "tbi1": vcf1 + ".tbi",
        "tbi2": vcf2 + ".tbi",
        "hmap": config["haplotype_maps"][m],
    }


rule pair_lod_block:
    input:
        unpack(pair_lod_inputs),
    output:
        parquet="results/pair_lod/{job}/{left}~{right}.parquet",
    conda:
        "../envs/snp_prioritization.yaml"
    log:
        "logs/pair_lod/{job}/{left}~{right}.log",
    script:
        "../scripts/pair_lod.py"


def flagged_pair_parquets(wildcards):
    # Force `cases` to run, then enumerate one parquet per flagged (UNEXPECTED_*) pair.
    tsv = checkpoints.cases.get().output.tsv
    df = pl.read_csv(tsv, separator="\t").filter(
        pl.col("RESULT").str.starts_with("UNEXPECTED_")
    )
    return sorted(
        {
            f"results/pair_lod/{r['COMPARISON']}/{r['LEFT_GROUP_VALUE']}~{r['RIGHT_GROUP_VALUE']}.parquet"
            for r in df.iter_rows(named=True)
        }
    )


rule pair_lod_all:
    input:
        flagged_pair_parquets,
    output:
        touch("results/pair_lod/.done"),
