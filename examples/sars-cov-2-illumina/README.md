# SARS-CoV-2 Illumina project

This is a real project configuration, not bundled sequence data.

Download a versioned Nextclade dataset before running:

```bash
nextclade dataset get --name sars-cov-2 --output-dir datasets/sars-cov-2
```

Record the dataset version from `datasets/sars-cov-2/pathogen.json`. Add paired FASTQ files,
edit `samples.csv`, validate, inspect the command plan, and run:

```bash
viroflow validate --config config.yaml
viroflow run --config config.yaml --dry-run
viroflow run --config config.yaml
```

For amplicon data, configure the exact primer BED used in the laboratory. Do not mix primer
schemes or reference coordinate systems.
