# Concatenate the per-job crosscheck TSVs into one parquet and emit the
# cross-participant pairs (different individuals) ranked by LOD.


rule aggregate:
    input:
        tsvs=expand(
            "results/crosscheck/crosscheck_{job}.tsv",
            job=list(config["crosscheck_jobs"].keys()),
        ),
    output:
        combined="results/combined.parquet",
        cross="results/cross_participant.tsv",
    conda:
        "../envs/polars.yaml"
    log:
        "logs/aggregate.log",
    script:
        "../scripts/aggregate.py"
