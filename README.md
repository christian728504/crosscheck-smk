# crosscheck-smk

Picard-fingerprinting **sample-swap / relatedness** analysis across four modalities
(WGS, ATAC-seq, WGBS, RNA-seq), WGS-anchored by default, as a reproducible Snakemake workflow.

For every sample it extracts a Picard fingerprint, then crosschecks each **WGS**
sample against every **other-modality** sample (`WGS × other`) to flag pairs whose
genotypes don't match their expected participant — i.e. swaps and unexpected
relatives.

## Pipeline

```
manifest (pre-flight)            results/manifest.parquet
   └─ fingerprint (per sample,    ExtractFingerprint → reheader → bgzip → tabix
      per map)                    -> results/fingerprints/{map}/{modality}/{id}.vcf.gz
        └─ crosscheck            Picard CrosscheckFingerprints, LEFT(INPUT) × RIGHT(SECOND_INPUT)
             └─ aggregate        results/combined.parquet, results/cross_participant.tsv
                  └─ cases       results/cases.tsv + results/cases/case_{n}.html  (checkpoint)
                       └─ pair_lod   results/pair_lod/{comparison}/{left}~{right}.parquet
                                     (per-block LOD Δ for each flagged pair)
```

## Layout

```
config/config.yaml              # samples come from the manifest; tunables live here
workflow/
  Snakefile                     # includes + rule all
  rules/{common,fingerprint,crosscheck,aggregate,cases,pair_lod}.smk
  envs/{picard,bcftools,polars,viz,snp_prioritization}.yaml   # per-rule conda envs
  scripts/{construct_manifest,write_sample_map,gather_crosscheck,aggregate,build_cases,pair_lod}.py
profiles/
  default/config.yaml           # local execution (also the auto-applied workflow profile)
  slurm/config.yaml             # SLURM execution
resources/                      # reference fasta(+.fai) + 3 haplotype maps (copied in)
results/                        # all outputs (manifest, fingerprints/, crosscheck/, *.parquet, cases/)
```

## Prerequisites

- **Runner env (`snakemake-base`).** Snakemake + the pre-flight script run from this
  env (it has `snakemake`, the SLURM executor plugin, `polars`, `pyyaml`). Create it
  once from the pinned spec at the repo root:
  ```bash
  mamba env create -f runtime-environment.yaml      # or: conda env create -f runtime-environment.yaml
  # update an existing one after the spec changes:
  mamba env update -n snakemake-base -f runtime-environment.yaml --prune
  ```
  (`runtime-environment.yaml` is `conda env export --no-builds` of `snakemake-base`;
  `name: snakemake-base` is baked in, so the env name matches what the commands below
  use.) Activate with `conda activate snakemake-base`, or prefix commands with
  `conda run -n snakemake-base ...`.
- **Per-rule tools** (picard, bcftools/htslib, polars, bokeh/graphviz) come from the
  `workflow/envs/*.yaml` conda envs, built automatically into `.conda/` on first run.
- **Resources (not in git — large).** The reference fasta (+`.fai`/`.dict`) and the
  three haplotype maps are too large to version. Fetch them into `resources/` once,
  from the repo root:
  ```bash
  curl -fsSL https://users.wenglab.org/ramirezc/share/crosscheck-smk/resources.tar.zst | tar --zstd -xf -
  ```
  (Picard itself comes from its conda env; the aligned BAM/CRAM inputs are read in
  place from `/data/projects/mohd/...` and are not bundled.)

## Usage

```bash
cd /zata/zippy/ramirezc/Projects/crosscheck-smk
conda activate snakemake-base      # or prefix commands with: conda run -n snakemake-base

# 0. One-time setup: fetch the large resource files (reference + haplotype maps) into
#    resources/ (not versioned in git). Run from the repo root.
curl -fsSL https://users.wenglab.org/ramirezc/share/crosscheck-smk/resources.tar.zst | tar --zstd -xf -

# 1. Pre-flight: build the manifest (samples + source paths) — run once, or whenever
#    the GCP-Upload data model changes.
python workflow/scripts/construct_manifest.py --config config/config.yaml --out results/manifest.parquet

# 2. Dry-run (always check the DAG first)
snakemake --profile profiles/default -n

# 3a. Run locally
snakemake --profile profiles/default

# 3b. Run on SLURM
snakemake --profile profiles/slurm
```

`config/config.yaml` is loaded automatically (declared via `configfile:` in the
workflow) — only pass `--configfile` to point at a *different* config. First run also
solves the conda envs; pre-build them with `--conda-create-envs-only`.

## Configuration (`config/config.yaml`)

| Key | Meaning |
|---|---|
| `crosscheck.scatter_chunk` | `null` = one Picard run per job over all WGS; `N` = split the WGS anchor into chunks of `N`, run each as its own job, then concatenate. **Required on SLURM** (see below). |
| `crosscheck.chunk_threads` / `chunk_gb` | Picard `--NUM_THREADS` / Java heap per chunk (scatter mode). |
| `crosscheck_threads` / `heaps.crosscheck_gb` | Threads / heap for the single un-scattered run (`scatter_chunk: null`). |
| `heaps.extract_gb` | Java heap for fingerprint extraction. |
| `lod_threshold` | LOD below which a pair is called a mismatch. |
| `crosscheck_jobs` | Pairwise comparisons + their haplotype map. `modalities: [A, B]` = A×B (cross); `[A]` = self (A×A). `crosscheck_map` is used for **both** extraction and crosscheck, so a sample is fingerprinted once per distinct map it appears in (e.g. WGS under full/filtered/coding_exons) — this maximizes the LOD blocks kept (no representative-SNP dropout). |

## Outputs (`results/`)

- `combined.parquet` — every configured crosscheck pair (all jobs, with a `COMPARISON` column).
- `cross_participant.tsv` — cross-participant pairs ranked by LOD (candidate relatives).
- `cases.tsv` — anomalous (`UNEXPECTED_*`) pairs plus the expected relationships they
  conflict with, annotated with `SUBGRAPH_ID`.
- `cases/case_{n}.html` — one interactive Bokeh graph per connected component
  (nodes labelled by participant, hover shows the MOHD accession).
- `pair_lod/{comparison}/{left}~{right}.parquet` — per-haplotype-block LOD
  contributions (`chrom, pos, name, maf, delta`) for each flagged (`UNEXPECTED_*`)
  pair in `cases.tsv`, via `snp-prioritization`'s `pair_lod` (`Σ delta == LOD_SCORE`).
  Fanned out by a checkpoint on `cases` (the pair set is unknown until `cases.tsv`
  exists); rank by `delta` to find the blocks driving a swap call.

## Notes

- **SLURM:** set `crosscheck.scatter_chunk` (e.g. 32) — the un-scattered run wants 128
  threads, which exceeds the 88-CPU node limit. The SLURM profile activates per-rule
  conda envs on compute nodes (`conda-base-path`) and disables the local-cores check
  (`core_check=false`). Keep `chunk_count × chunk_threads ≤ --cores` for parallel waves.
- **Picard exit code 1** (a flagged pair) is treated as success; only `>1` is a failure.
- Comparisons are **configurable** via `crosscheck_jobs` (any pairwise A×B, or self
  A×A). The default config is WGS-anchored. The `cases` diagnostic stays WGS-anchored
  (it only renders WGS-on-left cross-modality flags).

## Troubleshooting

First triage: the per-rule log (`logs/<rule>/.../<sample>.log`) has the tool's stderr;
the SLURM wrapper log (`.snakemake/slurm_logs/rule_<rule>/.../<jobid>.log`) has the
`host:` line and scheduler-side errors; `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS,NodeList`
shows what SLURM thinks happened.

### Jobs submit but never start (stuck `PENDING`)
- `squeue --me -o "%.18i %.8T %.20R"` → reason in the last column.
- **`(Priority)` / `(Resources)`** = the scheduler hasn't allocated yet (busy cluster,
  fair-share, or your jobs require scarce nodes). Check `squeue --me --start` for ETAs
  and `sinfo -p <partition> -o "%.14f %.6t %.5D"` for node availability.
- **A `constraint:` is starving you.** If a rule pins `constraint: cascadelake` but
  those nodes are down/drained while idle `broadwell` nodes sit free, jobs queue
  forever. Drop the constraint for light rules (they don't need a specific arch) — the
  data fits broadwell (88 CPU/500 GB).
- **`(JobHeldAdmin)`** = an admin/policy held your jobs (often triggered by submitting
  *thousands* of tiny jobs at once). They won't run until released; cancel + resubmit
  with fewer concurrent jobs (`jobs:` lower), throttle (`max-jobs-per-second`,
  `max-status-checks-per-second`), and ask the admin. Held jobs from a dead run are
  orphaned — `scancel` them rather than releasing.

### `MissingOutputException` / "completed successfully, but output files are missing"
- Almost always a **bad compute node missing the `/data` mount** (the input can't be
  read, the job exits 0 doing nothing). Check the job's SLURM log for
  `parent dir not present`, find the node from its `host:` line, and **exclude it**:
  `--exclude=zXXXX` in `profiles/slurm` `default-resources.slurm_extra`. Enumerate all
  bad nodes:
  ```bash
  grep -rl "parent dir not present" .snakemake/slurm_logs/ \
    | while read f; do grep -m1 '^host:' "$f"; done | awk '{print $2}' | sort | uniq -c
  ```
  Full write-up + node-level fix: **`docs/incident-data-mount.md`**.
- If it's genuine shared-FS latency (rare; output appears seconds later), raise
  `latency-wait`.

### One failure aborts the whole run
- Set **`keep-going: true`** so unrelated jobs finish and only the failed branch is
  left incomplete. `retries: N` re-runs transient failures — but note a deterministic
  failure (corrupt input, bad node that keeps grabbing the resubmit) will exhaust
  retries; fix the cause or exclude the node/sample.

### A specific sample always fails (`Invalid GZIP header` / `Invalid BGZF header`)
- The **input BAM/CRAM is corrupt** (a bad BGZF block; `samtools quickcheck` passes
  because it only checks header+EOF). Confirm with a region read:
  `samtools view -c <bam> chr1:1-1000000`. Fix the file upstream, and meanwhile drop it
  via `exclude_samples:` in `config.yaml` then rebuild the manifest.

### Picard `A reference dictionary is required for creating Tribble indices`
- `ExtractFingerprint` writes the VCF with htsjdk index-on-the-fly (independent of
  `--CREATE_INDEX`), which needs the reference **`.dict`** next to the fasta. Ensure
  `resources/<ref>.dict` exists (alongside `.fasta` + `.fai`).

### Snakemake won't build the DAG / parse errors
- `Manifest not found` → run the pre-flight `construct_manifest.py` first.
- Parse-time `ModuleNotFoundError: polars` → the **runner env** (`snakemake-base`)
  needs `polars` (the Snakefile reads the manifest at parse time).
- `--configfile` is **not required** — `config/config.yaml` is loaded via the
  `configfile:` directive; pass `--configfile` only to point at a different config.
- After a killed run, clear the lock with `snakemake --unlock` before resubmitting.

### Conda envs on compute nodes
- The SLURM profile sets `conda-base-path` (shared Miniforge3) so workers can activate
  per-rule envs, and `conda-prefix: .conda` so envs build once on the shared FS. No
  `--precommand` is needed (that flag is only for sourcing conda on remote workers).

### Harmless noise
- An `OrderedDict({...}) {...}` dump on every `--profile` run is a stray
  `print()` in the installed Snakemake's `profiles.py` (`ProfileConfigFileParser`),
  not this workflow. Ignore it (or comment that line in the env's `snakemake/profiles.py`).
