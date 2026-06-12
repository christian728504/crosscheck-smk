# Per-sample fingerprinting, split into two rules so each uses a minimal conda env:
#   extract_fingerprint -> Picard/Java produces a raw single-sample VCF (temp)
#   reheader_index      -> bcftools/htslib rename the sample to its dataset_id,
#                          then bgzip + tabix into the final .vcf.gz (+ .tbi)


rule extract_fingerprint:
    input:
        bam=extract_input_bam,
        ref=REF_FASTA,
        fai=REF_FASTA + ".fai",
        dict=REF_DICT,
        hmap=extract_map,
    output:
        vcf=temp("results/fingerprints/{map}/{modality}/{id}.vcf"),
        # Picard writes this Tribble index on the fly (htsjdk default); unused
        # downstream (we bgzip + tabix), so temp() reclaims it.
        idx=temp("results/fingerprints/{map}/{modality}/{id}.vcf.idx"),
    params:
        heap=config["heaps"]["extract_gb"],
    threads: 1
    resources:
        mem_mb=lambda wc: (config["heaps"]["extract_gb"] + 2) * 1024,
        runtime=120,
    conda:
        "../envs/picard.yaml"
    log:
        "logs/extract/{map}/{modality}/{id}.log",
    shell:
        r"""
        picard -Xmx{params.heap}g ExtractFingerprint \
            --INPUT {input.bam} \
            --OUTPUT {output.vcf} \
            --HAPLOTYPE_MAP {input.hmap} \
            --REFERENCE_SEQUENCE {input.ref} \
            --VALIDATION_STRINGENCY SILENT \
            > {log} 2>&1
        """


rule reheader_index:
    input:
        vcf="results/fingerprints/{map}/{modality}/{id}.vcf",
    output:
        gz="results/fingerprints/{map}/{modality}/{id}.vcf.gz",
        tbi="results/fingerprints/{map}/{modality}/{id}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=30,
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/reheader/{map}/{modality}/{id}.log",
    shell:
        r"""
        {{
          tmp="{input.vcf}.reheader.tmp"
          printf '%s\n' "{wildcards.id}" | bcftools reheader -s - "{input.vcf}" -o "$tmp"
          mv -f "$tmp" "{input.vcf}"
          bgzip -f "{input.vcf}"
          tabix -p vcf "{output.gz}"
        }} > {log} 2>&1
        """
