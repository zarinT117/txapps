# Scientific Scope

ViroFlow is designed for reproducible viral-genome surveillance and research screening.

## Appropriate uses

- raw-read processing to assembly and consensus
- haploid consensus variant calling for viral genomes
- primer-aware amplicon processing when the exact primer BED is configured
- minor-variant tabulation through iVar when sufficient depth is available
- configurable antigenic-site and escape-marker evidence summaries
- segmented-virus reassortment candidate screening
- cohort-level ML prioritization using validated labels

## Required researcher inputs

- virus-appropriate reference FASTA
- correct segment names and coding coordinates
- curated lineage references
- evidence-cited antigenic sites and escape markers
- local, versioned Nextclade datasets when Nextclade is used
- validation labels for supervised ML

## Interpretation boundaries

Sequence evidence alone does not prove antigenic drift, antigenic shift, vaccine escape,
transmissibility, virulence, severity, or vaccine failure. ViroFlow outputs should be interpreted
with laboratory evidence, phylogenetics, epidemiology, metadata review, and contamination checks.

## Publication checklist

- archive the ViroFlow version, commit hash, and `run_manifest.json`
- archive the analysis profile and all references used
- report sequencing platform, library method, primer scheme, and thresholds
- inspect coverage breadth and ambiguous bases before phenotype interpretation
- validate ML models on independent data from the intended use setting
