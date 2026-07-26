# ViroFlow

ViroFlow is an installable command-line workflow for viral whole-genome sequencing. It
connects raw-read QC and assembly to haploid variant calling, masked consensus generation,
genotyping, and configurable antigenic evidence screening.

The repository is pathogen-agnostic: the executable workflow is reusable, while biological
interpretation comes from a versioned **analysis profile** containing the correct reference,
coding coordinates, lineage references, antigenic sites, vaccine sequence, and evidence-linked
mutations for the virus being studied.

## What it produces

For each sample, `viroflow run` creates:

- QC-filtered reads and an HTML/JSON QC report
- de novo contigs (`MEGAHIT` for Illumina; `Flye` for Nanopore)
- reference-aligned BAM and all-position depth table
- raw and depth/quality-filtered haploid VCF
- a reference-guided consensus with low-depth regions masked as `N`
- optional Nextclade clade/lineage and QC output from a pinned local dataset
- JSON, TSV, and standalone HTML comparative-analysis reports
- a machine-readable run manifest containing every command and threshold

The comparative report includes nucleotide and amino-acid changes, configured antigenic-site
changes, best-reference segment genotype calls, a segmented-virus reassortment/shift screen,
and a transparent vaccine-escape prioritization score.

## Install after `git clone`

The full workflow runs on Linux or WSL with Conda/Mamba:

```bash
git clone https://github.com/zarinT117/txapps.git
cd txapps
mamba env create -f environment.yml
conda activate viroflow
viroflow --version
viroflow doctor
```

The sequence-analysis/reporting layer also works without the external bioinformatics tools:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

## Quick demonstration

The repository includes a tiny synthetic, two-segment fixture:

```bash
viroflow analyze \
  --profile examples/demo/profile.yaml \
  --input examples/demo/query.fasta \
  --sample synthetic-demo \
  --output demo-results
```

Open `demo-results/report.html`. The fixture deliberately gives one segment its closest match
in `lineage_A` and another in `lineage_B`, exercising the reassortment-candidate screen. It is
software test data, not a biological reference.

## Configure a real project

Create starter files:

```bash
viroflow init my-virus-project
cd my-virus-project
```

Edit `samples.csv`:

```csv
sample,r1,r2,platform
case01,reads/case01_R1.fastq.gz,reads/case01_R2.fastq.gz,illumina
case02,reads/case02.fastq.gz,,nanopore
```

Place a curated reference in `references/reference.fasta`, then edit `config.yaml` and
`profile.yaml`. Inspect the complete execution plan before spending compute:

```bash
viroflow run --config config.yaml --dry-run
viroflow run --config config.yaml
```

Paths in YAML and CSV files are resolved relative to the YAML file. FASTA record IDs must match
the profile's segment names. Coding coordinates and amino-acid positions are 1-based and
inclusive.

### Genotype references

`lineage_references` is a multi-record FASTA whose IDs use `lineage|segment`:

```text
>lineage_A|HA
...
>lineage_A|NA
...
>lineage_B|HA
...
```

ViroFlow reports the closest reference for each segment and its identity margin over the
runner-up. For segmented viruses, confident calls to multiple lineages are flagged as a
**candidate reassortment signal**. Confirm candidates with segment-specific phylogenetics,
read-level mixture/contamination checks, and epidemiologic context.

For supported pathogens, configure `workflow.nextclade_dataset` to a downloaded, versioned
Nextclade dataset. Pinning the local dataset makes genotype results reproducible:

```bash
nextclade dataset get --name '<dataset-name>' --output-dir datasets/my-dataset
```

### Antigenic and vaccine-escape profile

Profiles define antigenic positions and evidence-linked markers rather than baking a
single-virus marker list into the software:

```yaml
antigenic_sites:
  surface_protein: [123, 145, 156]

escape_markers:
  surface_protein:
    - mutation: A123T
      weight: 1.0
      evidence: "PMID:00000000"
```

Every marker must include a citation. Keep profiles under version control and review them when
surveillance or experimental evidence changes.

The reported escape score is an auditable triage index:

```text
70 × matched-marker-weight fraction
+ 20 × configured-antigenic-site change fraction
+ 10 × min(1, vaccine nucleotide distance / 0.05)
```

It is a ranking aid, not a probability, neutralization titre, vaccine-effectiveness estimate,
or clinical conclusion.

## Workflow outline

```text
FASTQ
  ├─ fastp (Illumina) / Filtlong (Nanopore)
  ├─ MEGAHIT / Flye ───────────────> de novo contigs
  └─ minimap2 → samtools → bcftools
                     ├──────────────> BAM + depth + VCF
                     └──────────────> low-depth-masked consensus
                                           ├─ Nextclade genotype/QC (optional)
                                           └─ ViroFlow comparative report
```

The commands follow the upstream interfaces documented by
[fastp](https://github.com/OpenGene/fastp),
[MEGAHIT](https://github.com/voutcn/megahit),
[minimap2](https://github.com/lh3/minimap2),
[SAMtools](https://www.htslib.org/doc/samtools.html),
[BCFtools](https://samtools.github.io/bcftools/bcftools), and
[Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/usage.html).

## Quality and interpretation

- Use a reference appropriate to the virus, segment, assay, and sampling period.
- Confirm primer-derived workflows have had primer sequences removed with an assay-aware method.
- Inspect coverage, allele balance, contamination, and ambiguous calls before interpreting a
  consensus.
- De novo contigs and reference-guided consensus answer different questions; the workflow keeps
  both.
- Sequence divergence alone does not establish antigenic drift, antigenic shift, immune escape,
  phenotype, transmissibility, severity, or vaccine effectiveness.
- Never use the synthetic demo profile for research samples.

## Development

```bash
python -m pip install -e ".[dev]"
ruff check src tests
pytest
```

Licensed under the [MIT License](LICENSE).
