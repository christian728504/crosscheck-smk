# Stage E: sample-swap cases. From the aggregated crosscheck results, build a table
# of anomalous (UNEXPECTED_*) comparisons plus the expected relationships they
# conflict with, and one interactive Bokeh graph per connected component.


rule cases:
    input:
        combined="results/combined.parquet",
        manifest=config["manifest"],
    output:
        tsv="results/cases.tsv",
        plots=directory("results/cases"),
    conda:
        "../envs/viz.yaml"
    log:
        "logs/cases.log",
    script:
        "../scripts/build_cases.py"
