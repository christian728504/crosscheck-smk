# Shared setup: load config + manifest at parse time, derive sample lookups, and
# declare wildcard constraints used across the workflow.
#
# The manifest is produced out-of-band by workflow/scripts/construct_manifest.py,
# so the full set of samples is known before snakemake builds the DAG -- a plain
# parse-time read is correct here and a checkpoint would only add a redundant
# DAG re-evaluation.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger


configfile: "config/config.yaml"


_MANIFEST = Path(config["manifest"])
if not _MANIFEST.exists():
    raise WorkflowError(
        f"Manifest not found: {_MANIFEST}\n"
        "Build it first with the pre-flight script, e.g.:\n"
        "  conda run -n snakemake-base python workflow/scripts/construct_manifest.py "
        f"--config config/config.yaml --out {_MANIFEST}"
    )

# Belt-and-suspenders: the script already keeps only existing files.
manifest = pl.read_parquet(_MANIFEST).filter(pl.col("file_exists"))

# data_modality label -> output subdir (wgs/atac/wgbs/rna).
SUBDIR = {label: spec["subdir"] for label, spec in config["modalities"].items()}

# data_modality label -> list of dataset_ids present in the manifest.
IDS = {
    label: manifest.filter(pl.col("data_modality") == label)
    .get_column("entity:dataset_id")
    .to_list()
    for label in config["modalities"]
}

# dataset_id -> source alignment path (each dataset_id belongs to one modality).
FILE_OF = {
    row["entity:dataset_id"]: row["file"] for row in manifest.iter_rows(named=True)
}

# dataset_id -> participant_id, for the SAMPLE_INDIVIDUAL_MAP.
PARTICIPANT_OF = dict(
    zip(
        manifest.get_column("entity:dataset_id").to_list(),
        manifest.get_column("participant_id").to_list(),
    )
)

# Fail loudly rather than run a degenerate one-modality crosscheck.
for _job, _spec in config["crosscheck_jobs"].items():
    _empty = [m for m in _spec["modalities"] if not IDS.get(m)]
    if _empty:
        raise WorkflowError(
            f"Crosscheck job '{_job}' has modalities with zero samples: {_empty}. "
            "Check the manifest or drop the job from config."
        )


REF_FASTA = config["reference"]["fasta"]
# Sequence dictionary, auto-discovered by Picard next to the fasta. ExtractFingerprint
# writes the VCF with htsjdk index-on-the-fly (independent of --CREATE_INDEX), which
# needs this .dict -- without it Picard fails ("reference dictionary is required ...").
REF_DICT = REF_FASTA.removesuffix(".fasta") + ".dict"


# --- Crosscheck scatter-gather (pairwise, Picard SECOND_INPUT) -------------------
# Each job compares a LEFT modality (Picard --INPUT) against a RIGHT modality
# (--SECOND_INPUT). modalities: [A, B] -> A x B (cross); [A] -> A x A (self). The
# LEFT modality is scattered into chunks of `scatter_chunk` samples; each chunk runs
# against ALL of RIGHT. Chunk counts are deterministic at parse time.
_cc = config["crosscheck"]
SCATTER = int(_cc["scatter_chunk"])
if SCATTER < 1:
    raise WorkflowError(f"crosscheck.scatter_chunk must be a positive integer, got {SCATTER}")

# dataset_id -> modality, for building the per-sample fingerprint path.
MODALITY_OF = dict(
    zip(
        manifest.get_column("entity:dataset_id").to_list(),
        manifest.get_column("data_modality").to_list(),
    )
)


def job_pair(job):
    """(left_modality, right_modality) for a job; right == left for a self comparison."""
    mods = config["crosscheck_jobs"][job]["modalities"]
    if len(mods) not in (1, 2):
        raise WorkflowError(
            f"Crosscheck job '{job}': modalities must have 1 (self) or 2 entries, got {mods}."
        )
    return mods[0], (mods[1] if len(mods) == 2 else mods[0])


def is_self(job):
    left, right = job_pair(job)
    return left == right


def left_chunks(job):
    """LEFT modality dataset_ids split into scatter chunks of `scatter_chunk` each."""
    ids = IDS[job_pair(job)[0]]
    return [ids[i : i + SCATTER] for i in range(0, len(ids), SCATTER)] or [[]]


def n_chunks(job):
    return len(left_chunks(job))


def fp(dataset_id, map_name):
    """Reheadered fingerprint VCF path for a sample extracted with a given map.
    Keyed by (map, modality) so a sample used in comparisons with different maps is
    extracted once per map (extraction map == that comparison's crosscheck map)."""
    return f"results/fingerprints/{map_name}/{SUBDIR[MODALITY_OF[dataset_id]]}/{dataset_id}.vcf.gz"


def job_map(job):
    """Haplotype map name for a comparison -- used for BOTH extraction and crosscheck."""
    return config["crosscheck_jobs"][job]["crosscheck_map"]


# Advisory check (local execution only): the largest job's single wave of chunks
# should fit in the available cores so they run in parallel. Snakemake still caps
# concurrency by cores, so this only warns. Disabled under SLURM via core_check=false.
_core_check = str(config.get("core_check", True)).strip().lower() not in {"false", "0", "no", "off", "none"}
if _core_check:
    _chunk_threads = int(_cc["chunk_threads"])
    _max_chunks = max((n_chunks(j) for j in config["crosscheck_jobs"]), default=0)
    try:
        _cores = int(workflow.cores)
    except (TypeError, ValueError, NameError):
        _cores = None
    if _cores and _max_chunks * _chunk_threads > _cores:
        logger.warning(
            f"crosscheck scatter: up to {_max_chunks} chunks x {_chunk_threads} threads "
            f"= {_max_chunks * _chunk_threads} > {_cores} cores; chunks will run in waves "
            f"(reduce crosscheck.chunk_threads or raise scatter_chunk to parallelize fully)."
        )


wildcard_constraints:
    # dataset_ids contain underscores/dots/dashes but never a path separator.
    id=r"[^/]+",
    map="|".join(config["haplotype_maps"].keys()),
    modality="|".join(SUBDIR.values()),
    job="|".join(config["crosscheck_jobs"].keys()),
    chunk=r"\d+",


def extract_input_bam(wildcards):
    """Source BAM/CRAM for one dataset_id."""
    return FILE_OF[wildcards.id]


def extract_map(wildcards):
    """Haplotype map path for an extraction, from the {map} wildcard (== the
    crosscheck map of the comparison this fingerprint feeds)."""
    return config["haplotype_maps"][wildcards.map]


# Every (map, dataset_id) fingerprint a comparison needs: each job extracts both of
# its modalities with the job's map. Deduped, so a sample shared by same-map jobs is
# extracted once.
FP_TARGETS = {
    (job_map(_job), _i)
    for _job in config["crosscheck_jobs"]
    for _mod in set(job_pair(_job))
    for _i in IDS[_mod]
}


def all_fingerprints():
    """Every reheadered, bgzipped fingerprint VCF the comparisons imply."""
    return [fp(i, m) for (m, i) in FP_TARGETS]
