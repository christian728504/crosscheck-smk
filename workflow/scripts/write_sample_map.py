"""Write a Picard SAMPLE_INDIVIDUAL_MAP for one crosscheck job.

Two-column TSV (no header): dataset_id <TAB> participant_id, one line per input
sample. Same-participant pairs become EXPECTED_MATCH in CrosscheckFingerprints.
Run via Snakemake's ``script:`` directive, so the ``snakemake`` object is injected.
"""

from pathlib import Path

ids = snakemake.params.ids  # noqa: F821 - injected by Snakemake
participant_of = snakemake.params.participant_of  # noqa: F821

seen = set()
lines = []
for dataset_id in ids:
    if dataset_id in seen:
        continue
    seen.add(dataset_id)
    lines.append(f"{dataset_id}\t{participant_of[dataset_id]}")

text = "\n".join(lines)
Path(snakemake.output.tsv).write_text(text + "\n" if text else "")  # noqa: F821
