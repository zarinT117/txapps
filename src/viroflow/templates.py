CONFIG_TEMPLATE = """\
workflow:
  samples: samples.csv
  reference_fasta: references/reference.fasta
  output_dir: results
  platform: illumina
  threads: 4
  min_read_length: 50
  min_depth: 10
  min_variant_quality: 20
  # Download and pin a Nextclade dataset, then set its local directory here.
  # nextclade_dataset: datasets/nextclade
  # Optional production controls:
  # host_reference_fasta: references/host.fasta
  # primer_bed: references/primers.bed
  # annotation_gff: references/annotation.gff
  # primer_min_length: 50
  # primer_min_quality: 20

analysis:
  profile: profile.yaml

minor_variants:
  enabled: false
  min_frequency: 0.03
  min_depth: 100
  min_base_quality: 20
"""

SARS_COV_2_CONFIG_TEMPLATE = """\
workflow:
  samples: samples.csv
  reference_fasta: references/sars-cov-2.fasta
  output_dir: results
  platform: illumina
  threads: 8
  min_read_length: 50
  min_depth: 10
  min_variant_quality: 20
  # Download/pin with:
  # nextclade dataset get --name sars-cov-2 --output-dir datasets/nextclade/sars-cov-2
  nextclade_dataset: datasets/nextclade/sars-cov-2
  # For ARTIC or tiled-amplicon data, set the matching primer BED.
  # primer_bed: references/artic_primers.bed
  # Optional host depletion against a licensed/local human reference.
  # host_reference_fasta: references/human_or_host.fasta

analysis:
  profile: profile.yaml

minor_variants:
  enabled: false
  min_frequency: 0.03
  min_depth: 100
  min_base_quality: 20
"""

SAMPLES_TEMPLATE = """\
sample,r1,r2,platform
sample01,reads/sample01_R1.fastq.gz,reads/sample01_R2.fastq.gz,illumina
"""

PROFILE_TEMPLATE = """\
# Segment names must match FASTA record IDs.
name: example-virus
segmented: false
reference_fasta: references/reference.fasta
# vaccine_reference_fasta: references/vaccine.fasta
# FASTA IDs use lineage|segment, for example: lineage_A|genome
# lineage_references: references/lineages.fasta
lineage_min_margin: 0.002

coding_regions:
  - name: surface_protein
    segment: genome
    start: 1
    end: 300
    strand: 1

# 1-based amino-acid positions in each protein.
antigenic_sites:
  surface_protein: []

# Every configured marker requires a traceable citation.
escape_markers:
  surface_protein: []
  # - mutation: A123T
  #   weight: 1.0
  #   evidence: "PMID, DOI, or stable public-health source"
"""
