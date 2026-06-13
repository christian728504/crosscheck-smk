# Per-block LOD decomposition for the sample pairs in cases.tsv (checkpoint fan-out).
#
# `cases` is a checkpoint: the set of pairs it keeps is only known once cases.tsv
# exists, so the number of per-pair parquets is dynamic. For every row in cases.tsv
# (UNEXPECTED_* flags and the EXPECTED_MATCH rows they conflict with) we run
# snp_prioritization.pair_lod on that pair's two fingerprint VCFs and the job's
# haplotype map, writing the per-block deltas to
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


def case_pair_parquets(wildcards):
    # Force `cases` to run, then enumerate one parquet per pair in cases.tsv -- every
    # comparison it kept (the UNEXPECTED_* flags AND the EXPECTED_MATCH rows they
    # conflict with), so the block-level evidence is available for both sides.
    tsv = checkpoints.cases.get().output.tsv
    df = pl.read_csv(tsv, separator="\t")
    return sorted(
        {
            f"results/pair_lod/{r['COMPARISON']}/{r['LEFT_GROUP_VALUE']}~{r['RIGHT_GROUP_VALUE']}.parquet"
            for r in df.iter_rows(named=True)
        }
    )


rule pair_lod_all:
    input:
        case_pair_parquets,
    output:
        touch("results/pair_lod/.done"),
