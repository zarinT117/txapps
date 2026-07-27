# ViroFlow

ViroFlow is an installable viral-genome analysis workflow by **Tasnim Zarin**. It connects raw
reads to read QC, host depletion, de novo assembly, alignment, variant calling, consensus
generation, genotype evidence, antigenic drift/shift screening, vaccine-escape prioritization,
and optional machine-learning triage.

The software is pathogen-agnostic. Biological interpretation comes from versioned analysis
profiles supplied by the researcher: reference sequences, coding coordinates, lineage references,
antigenic sites, vaccine references, and evidence-cited escape markers.

## Install after git clone

Full raw-read workflows require Linux or WSL with Conda/Mamba:

```bash
git clone https://github.com/zarinT117/txapps.git
cd txapps
mamba env create -f environment.yml
conda activate viroflow
viroflow --version
viroflow doctor
```

For analysis/reporting and ML utilities only:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[ml]"
```

Every normal command prints the runtime banner:

```text
ViroFlow 0.2.0 | Author: Tasnim Zarin <tasnim.2001040@bau.edu.bd>
```

## What it produces

For each sample, `viroflow run` writes:

- filtered reads and upstream QC files
- optional host-depleted reads
- MEGAHIT or Flye de novo contigs
- minimap2/samtools BAM, index, flagstat, and all-position depth
- BCFtools raw, normalized, filtered, indexed VCF plus VCF stats
- low-depth-masked consensus FASTA
- optional iVar primer-trimmed BAM and minor-variant TSV
- optional Nextclade output from a pinned local dataset
- JSON, TSV, and standalone HTML antigenic evidence reports
- QC summary, cohort summary, and a provenance manifest with command history, thresholds,
  input SHA256 hashes, tool paths, and tool versions

## Create a project

```bash
viroflow init my-virus-project
cd my-virus-project
```

For SARS-CoV-2 Illumina surveillance scaffolding:

```bash
viroflow init sars2-project --preset sars-cov-2
```

Edit `samples.csv`:

```csv
sample,r1,r2,platform
case01,reads/case01_R1.fastq.gz,reads/case01_R2.fastq.gz,illumina
case02,reads/case02.fastq.gz,,nanopore
```

Validate and run:

```bash
viroflow validate --config config.yaml
viroflow run --config config.yaml --dry-run
viroflow run --config config.yaml --resume
```

Paths in YAML and CSV files are resolved relative to the YAML file. FASTA record IDs must match
profile segment names. Coordinates and amino-acid positions are 1-based and inclusive.

## Genotyping, drift, shift, and escape evidence

`lineage_references` is a FASTA with IDs in `lineage|segment` form. ViroFlow assigns each segment
to the closest configured reference and reports the identity margin. Segmented viruses with
confident segment calls from multiple lineages are flagged as candidate reassortment signals that
need phylogenetic, read-level, contamination, and epidemiologic confirmation.

Profiles also define antigenic sites and evidence-linked escape markers:

```yaml
antigenic_sites:
  surface_protein: [123, 145, 156]

escape_markers:
  surface_protein:
    - mutation: A123T
      weight: 1.0
      evidence: "PMID:00000000"
```

The escape score is an auditable triage index:

```text
70 * matched-marker-weight fraction
+ 20 * configured-antigenic-site change fraction
+ 10 * min(1, vaccine nucleotide distance / 0.05)
```

It is not a neutralization estimate, vaccine-effectiveness model, clinical conclusion, or public
health decision rule.

## Machine learning

ViroFlow includes real ML utilities for cohort-level prioritization:

```bash
viroflow ml features --results results --output models/cohort_features.csv

viroflow ml train \
  --features labelled_features.csv \
  --label-column escape_label \
  --model models/escape_classifier.joblib \
  --report models/escape_classifier.metrics.json

viroflow ml predict \
  --model models/escape_classifier.joblib \
  --features models/cohort_features.csv \
  --output models/escape_predictions.csv

viroflow ml anomaly \
  --features models/cohort_features.csv \
  --output models/anomaly_scores.csv
```

The supervised classifier uses imputation, scaling, balanced logistic regression, stratified
cross-validation, and saved schema metadata. The anomaly mode uses Isolation Forest for unusual
genome/QC profiles and ignores common metadata columns such as `label`; pass repeated
`--exclude-column` values for other non-feature columns. Use externally validated labels for escape,
drift, phenotype, or vaccine breakthrough modelling; ViroFlow will not invent labels from sequence
names.

See [docs/MODEL_CARD.md](docs/MODEL_CARD.md) and [docs/SCIENTIFIC_SCOPE.md](docs/SCIENTIFIC_SCOPE.md).

## Reproducible datasets

For supported pathogens, download and pin a local Nextclade dataset, then set
`workflow.nextclade_dataset` to that directory:

```bash
nextclade dataset get --name sars-cov-2 --output-dir datasets/nextclade/sars-cov-2
```

Nextclade datasets are versioned; keep the downloaded dataset with your analysis record. ViroFlow
uses the official Nextclade CLI dataset/run interface, iVar primer/minor-variant workflow, and
BCFtools mpileup/call/consensus workflow.

Primary documentation:
[Nextclade datasets](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html),
[Nextclade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/usage.html),
[iVar manual](https://andersen-lab.github.io/ivar/html/manualpage.html), and
[BCFtools manual](https://samtools.github.io/bcftools/bcftools).

## Development

```bash
python -m pip install -e ".[dev]"
ruff check src tests
pytest
python -m build
twine check dist/*
```

Synthetic FASTA files live only under `tests/fixtures/synthetic` and are for software tests, not
research interpretation.

Licensed under the [MIT License](LICENSE). Cite this software with [CITATION.cff](CITATION.cff).
