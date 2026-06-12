# Pairwise crosschecks, scatter-gather over the LEFT modality.
#
# Each job compares a LEFT modality (Picard --INPUT) against a RIGHT modality
# (--SECOND_INPUT) with CHECK_ALL_OTHERS -> exactly LEFT x RIGHT pairs. The LEFT
# modality is split into chunks (left_chunks(job), from config["crosscheck"]
# ["scatter_chunk"]): each chunk runs its LEFT subset against ALL of RIGHT, and
# crosscheck_gather concatenates the per-chunk TSVs into crosscheck_{job}.tsv.
#
# Self comparison: modalities: [A] -> LEFT == RIGHT == A. The chunk's samples then
# appear in both --INPUT and --SECOND_INPUT, so Picard emits self-pairs and both
# orientations; the gather dedups (drops self-pairs, keeps one row per unordered pair).


def job_haplotype_map(wildcards):
    # Same map used to extract this job's fingerprints (job_map from common.smk).
    return config["haplotype_maps"][job_map(wildcards.job)]


def cc_threads():
    return int(config["crosscheck"]["chunk_threads"])


def cc_heap_gb():
    return int(config["crosscheck"]["chunk_gb"])


def chunk_left_vcfs(wildcards):
    m = job_map(wildcards.job)
    return [fp(i, m) for i in left_chunks(wildcards.job)[int(wildcards.chunk)]]


def chunk_right_vcfs(wildcards):
    m = job_map(wildcards.job)
    return [fp(i, m) for i in IDS[job_pair(wildcards.job)[1]]]


def chunk_map_ids(wildcards):
    # LEFT chunk + all of RIGHT (write_sample_map.py dedups, so self-overlap is fine).
    return left_chunks(wildcards.job)[int(wildcards.chunk)] + IDS[job_pair(wildcards.job)[1]]


def gather_inputs(wildcards):
    return [
        f"results/crosscheck/chunks/{wildcards.job}/chunk_{c}.tsv"
        for c in range(n_chunks(wildcards.job))
    ]


rule crosscheck_chunk_map:
    # SAMPLE_INDIVIDUAL_MAP for this chunk: dataset_id -> participant_id for every
    # sample in the run (the LEFT chunk + all of RIGHT).
    output:
        tsv="results/crosscheck/chunks/{job}/chunk_{chunk}.sample_map.tsv",
    params:
        ids=chunk_map_ids,
        participant_of=PARTICIPANT_OF,
    log:
        "logs/crosscheck/{job}/chunk_{chunk}.sample_map.log",
    script:
        "../scripts/write_sample_map.py"


rule crosscheck_chunk:
    input:
        left=chunk_left_vcfs,
        right=chunk_right_vcfs,
        left_tbi=lambda wc: [v + ".tbi" for v in chunk_left_vcfs(wc)],
        right_tbi=lambda wc: [v + ".tbi" for v in chunk_right_vcfs(wc)],
        smap="results/crosscheck/chunks/{job}/chunk_{chunk}.sample_map.tsv",
        ref=REF_FASTA,
        fai=REF_FASTA + ".fai",
        hmap=job_haplotype_map,
    output:
        tsv="results/crosscheck/chunks/{job}/chunk_{chunk}.tsv",
    params:
        inputs=lambda wc, input: " ".join(f"--INPUT {v}" for v in input.left),
        second=lambda wc, input: " ".join(f"--SECOND_INPUT {v}" for v in input.right),
        heap=cc_heap_gb(),
        lod=config["lod_threshold"],
    threads: cc_threads()
    resources:
        mem_mb=lambda wc: (cc_heap_gb() + 8) * 1024,
        runtime=1440,
    conda:
        "../envs/picard.yaml"
    log:
        "logs/crosscheck/{job}/chunk_{chunk}.log",
    shell:
        # INPUT=LEFT chunk, SECOND_INPUT=all RIGHT -> LEFT x RIGHT pairs.
        # Picard exits 1 when any pair is flagged (normal); only >1 is a real failure.
        r"""
        set +e
        picard -Xmx{params.heap}g CrosscheckFingerprints \
            {params.inputs} \
            {params.second} \
            --SAMPLE_INDIVIDUAL_MAP {input.smap} \
            --REFERENCE_SEQUENCE {input.ref} \
            --HAPLOTYPE_MAP {input.hmap} \
            --CROSSCHECK_BY SAMPLE \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --LOD_THRESHOLD {params.lod} \
            --CALCULATE_TUMOR_AWARE_RESULTS false \
            --NUM_THREADS {threads} \
            --OUTPUT {output.tsv} \
            > {log} 2>&1
        rc=$?
        if [ "$rc" -gt 1 ]; then
            echo "Picard CrosscheckFingerprints failed with exit code $rc" >> {log}
            exit "$rc"
        fi
        exit 0
        """


rule crosscheck_gather:
    # Concatenate the per-chunk TSVs into one per-job result. Cross jobs: pure concat
    # (each LEFT sample is in exactly one chunk). Self jobs: dedup (drop self-pairs,
    # keep one row per unordered pair) since INPUT/SECOND_INPUT overlap.
    input:
        chunks=gather_inputs,
    output:
        tsv="results/crosscheck/crosscheck_{job}.tsv",
    params:
        dedup=lambda wc: is_self(wc.job),
    conda:
        "../envs/polars.yaml"
    log:
        "logs/crosscheck/{job}.gather.log",
    script:
        "../scripts/gather_crosscheck.py"
