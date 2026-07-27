# ViroFlow ML Model Card

## Intended use

ViroFlow ML utilities support research triage of viral-genome cohorts after the standard
bioinformatics workflow has produced consensus, QC, antigenic evidence, and provenance outputs.

Supported modes:

- `viroflow ml features`: converts ViroFlow sample result directories into a numeric feature table
- `viroflow ml train`: trains a labelled classifier from researcher-provided labels
- `viroflow ml predict`: applies a saved classifier to matching feature tables
- `viroflow ml anomaly`: scores unusual cohort profiles with Isolation Forest

## Model design

The supervised model is a scikit-learn pipeline:

- `SimpleImputer`
- `StandardScaler`
- `LogisticRegression(class_weight="balanced", random_state=42)`

Metrics are estimated with stratified cross-validation. Saved model bundles include feature schema,
label column, classes, cross-validation metrics, ViroFlow version, author, and timestamp.

The anomaly model uses `IsolationForest(n_estimators=500, random_state=42)`.

## Inputs

Features include genome identity, N content, callable fraction, nucleotide and amino-acid change
counts, antigenic-site changes, matched marker weight, drift index, escape score, reassortment
screen flags, assembly metrics, and coverage metrics.

## Labels

Supervised models require external labels such as experimentally validated immune escape,
neutralization reduction class, vaccine-breakthrough status, or expert-reviewed drift priority.
Do not derive training labels from ViroFlow heuristic scores if the model is meant to validate the
same biological endpoint.

## Limitations

These models are research prioritization tools. They do not establish phenotype, antigenicity,
transmissibility, vaccine effectiveness, severity, or clinical outcome. Performance depends on the
virus, sampling frame, assay, label quality, class balance, and temporal validation.

## Recommended reporting

When publishing results, report the label definition, sample inclusion criteria, sequencing
platform, reference/profile versions, feature columns, cross-validation scheme, held-out temporal
or geographic validation when available, class balance, metrics, and model file hash.
