from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __author__, __version__


class MLError(ValueError):
    """Raised for user-correctable machine-learning input problems."""


def _require_sklearn() -> dict[str, Any]:
    try:
        import joblib
        import numpy as np
        from sklearn.ensemble import IsolationForest
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import (
            balanced_accuracy_score,
            f1_score,
            precision_score,
            recall_score,
            roc_auc_score,
        )
        from sklearn.model_selection import StratifiedKFold, cross_val_predict
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise MLError(
            "ML commands require optional dependencies. Install with "
            "`pip install -e .[ml]` or use the provided conda environment."
        ) from exc
    return {
        "IsolationForest": IsolationForest,
        "LogisticRegression": LogisticRegression,
        "Pipeline": Pipeline,
        "SimpleImputer": SimpleImputer,
        "StandardScaler": StandardScaler,
        "StratifiedKFold": StratifiedKFold,
        "balanced_accuracy_score": balanced_accuracy_score,
        "cross_val_predict": cross_val_predict,
        "f1_score": f1_score,
        "joblib": joblib,
        "np": np,
        "precision_score": precision_score,
        "recall_score": recall_score,
        "roc_auc_score": roc_auc_score,
    }


def _to_float(value: str, column: str, row_number: int) -> float:
    import math

    text = value.strip()
    if text == "":
        return math.nan
    try:
        return float(text)
    except ValueError as exc:
        raise MLError(f"non-numeric value in feature column {column!r} on row {row_number}") from exc


def read_feature_table(
    path: Path,
    label_column: str | None = None,
    sample_column: str = "sample",
    exclude_columns: set[str] | None = None,
) -> tuple[list[str], list[str], list[list[float]], list[str] | None]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise MLError(f"feature table has no header: {path}")
        fieldnames = list(reader.fieldnames)
        if sample_column not in fieldnames:
            raise MLError(f"feature table requires a {sample_column!r} column")
        if label_column and label_column not in fieldnames:
            raise MLError(f"label column {label_column!r} not found")
        reserved = {sample_column, label_column, *(exclude_columns or set())}
        feature_columns = [name for name in fieldnames if name not in reserved]
        if not feature_columns:
            raise MLError("feature table must contain at least one numeric feature column")

        sample_ids: list[str] = []
        features: list[list[float]] = []
        labels: list[str] = []
        for row_number, row in enumerate(reader, start=2):
            sample = (row.get(sample_column) or "").strip()
            if not sample:
                raise MLError(f"missing sample id on row {row_number}")
            sample_ids.append(sample)
            features.append([_to_float(row.get(column, ""), column, row_number) for column in feature_columns])
            if label_column:
                label = (row.get(label_column) or "").strip()
                if not label:
                    raise MLError(f"missing label in column {label_column!r} on row {row_number}")
                labels.append(label)
    return sample_ids, feature_columns, features, labels if label_column else None


def _flatten_sample_features(sample_dir: Path) -> dict[str, Any]:
    analysis_path = sample_dir / "05_report" / "analysis.json"
    qc_path = sample_dir / "05_report" / "qc_summary.json"
    if not analysis_path.exists() and not qc_path.exists():
        raise MLError(f"no ViroFlow report files found for sample directory: {sample_dir}")

    analysis = json.loads(analysis_path.read_text(encoding="utf-8")) if analysis_path.exists() else {}
    qc = json.loads(qc_path.read_text(encoding="utf-8")) if qc_path.exists() else {}
    segments = analysis.get("segments", {})
    proteins = analysis.get("proteins", {})
    segment_values = list(segments.values())

    identities = [float(item.get("identity", 0)) for item in segment_values]
    n_content = [float(item.get("n_content", 0)) for item in segment_values]
    callable_fraction = [
        float(item.get("callable_sites", 0)) / float(item.get("reference_length", 1) or 1)
        for item in segment_values
    ]
    matched_markers = [
        marker
        for protein in proteins.values()
        for marker in protein.get("matched_markers", [])
    ]
    antigenic_changes = [
        change
        for protein in proteins.values()
        for change in protein.get("antigenic_site_changes", [])
    ]
    aa_changes = [
        change
        for protein in proteins.values()
        for change in protein.get("amino_acid_changes", [])
    ]
    genotype_calls = analysis.get("genotype", {}).get("segment_calls", {})
    confident_lineages = {
        call.get("lineage")
        for call in genotype_calls.values()
        if call.get("confident") and call.get("lineage")
    }

    coverage = qc.get("coverage", {})
    assembly = qc.get("assembly", {})
    return {
        "sample": sample_dir.name,
        "segments": len(segment_values),
        "mean_identity": sum(identities) / len(identities) if identities else "",
        "min_identity": min(identities) if identities else "",
        "max_n_content": max(n_content) if n_content else "",
        "mean_callable_fraction": sum(callable_fraction) / len(callable_fraction)
        if callable_fraction
        else "",
        "nucleotide_changes_total": sum(len(item.get("nucleotide_changes", [])) for item in segment_values),
        "inserted_bases_total": sum(int(item.get("inserted_bases", 0)) for item in segment_values),
        "amino_acid_changes_total": len(aa_changes),
        "antigenic_site_changes_total": len(antigenic_changes),
        "matched_escape_markers_total": len(matched_markers),
        "matched_escape_marker_weight": sum(float(item.get("weight", 0)) for item in matched_markers),
        "drift_index": analysis.get("drift", {}).get("evidence_index", ""),
        "amino_acid_divergence": analysis.get("drift", {}).get("amino_acid_divergence", ""),
        "escape_score": analysis.get("vaccine_escape_screen", {}).get("priority_score", ""),
        "candidate_reassortment_signal": int(
            bool(analysis.get("shift_screen", {}).get("candidate_reassortment_signal", False))
        ),
        "confident_lineage_count": len(confident_lineages),
        "assembly_contigs": assembly.get("contigs", ""),
        "assembly_total_bases": assembly.get("total_bases", ""),
        "assembly_n50": assembly.get("n50", ""),
        "assembly_gc_fraction": assembly.get("gc_fraction", ""),
        "mean_depth": coverage.get("mean_depth", ""),
        "median_depth": coverage.get("median_depth", ""),
        "min_depth": coverage.get("min_depth", ""),
        "breadth_at_1x": coverage.get("breadth_at_1x", ""),
        "breadth_at_10x": coverage.get("breadth_at_10x", ""),
        "breadth_at_100x": coverage.get("breadth_at_100x", ""),
        "coverage_uniformity": coverage.get("uniformity_fraction_0_2x_to_2x_mean", ""),
    }


def extract_features(results_dir: Path, output_path: Path) -> Path:
    sample_dirs = sorted(path for path in results_dir.iterdir() if path.is_dir() and not path.name.startswith("_"))
    if not sample_dirs:
        raise MLError(f"no sample result directories found in {results_dir}")
    rows = [_flatten_sample_features(sample_dir) for sample_dir in sample_dirs]
    fieldnames = ["sample"] + sorted(name for name in rows[0] if name != "sample")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def train_classifier(
    features_path: Path,
    label_column: str,
    model_path: Path,
    report_path: Path,
    sample_column: str = "sample",
) -> dict[str, Any]:
    sk = _require_sklearn()
    np = sk["np"]
    sample_ids, feature_columns, features, labels = read_feature_table(features_path, label_column, sample_column)
    assert labels is not None
    classes, counts = np.unique(labels, return_counts=True)
    if len(classes) < 2:
        raise MLError("supervised training requires at least two label classes")
    min_class_count = int(counts.min())
    if min_class_count < 2:
        raise MLError("each label class needs at least two samples for cross-validation")

    pipeline = sk["Pipeline"](
        [
            ("imputer", sk["SimpleImputer"]()),
            ("scaler", sk["StandardScaler"]()),
            (
                "classifier",
                sk["LogisticRegression"](
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=42,
                ),
            ),
        ]
    )
    x_matrix = np.asarray(features, dtype=float)
    y_vector = np.asarray(labels)
    cv = sk["StratifiedKFold"](n_splits=min(5, min_class_count), shuffle=True, random_state=42)
    predicted = sk["cross_val_predict"](pipeline, x_matrix, y_vector, cv=cv)
    predicted_proba = sk["cross_val_predict"](pipeline, x_matrix, y_vector, cv=cv, method="predict_proba")
    metrics = {
        "samples": len(sample_ids),
        "features": len(feature_columns),
        "classes": {str(label): int(count) for label, count in zip(classes, counts)},
        "cross_validation_folds": cv.get_n_splits(),
        "balanced_accuracy": round(float(sk["balanced_accuracy_score"](y_vector, predicted)), 6),
        "macro_f1": round(float(sk["f1_score"](y_vector, predicted, average="macro", zero_division=0)), 6),
        "macro_precision": round(
            float(sk["precision_score"](y_vector, predicted, average="macro", zero_division=0)), 6
        ),
        "macro_recall": round(float(sk["recall_score"](y_vector, predicted, average="macro", zero_division=0)), 6),
    }
    try:
        if len(classes) == 2:
            metrics["roc_auc"] = round(float(sk["roc_auc_score"](y_vector, predicted_proba[:, 1])), 6)
        else:
            metrics["roc_auc_ovr"] = round(
                float(sk["roc_auc_score"](y_vector, predicted_proba, multi_class="ovr")), 6
            )
    except ValueError:
        metrics["roc_auc"] = None

    pipeline.fit(x_matrix, y_vector)
    bundle = {
        "model": pipeline,
        "feature_columns": feature_columns,
        "label_column": label_column,
        "sample_column": sample_column,
        "classes": [str(item) for item in classes],
        "metrics": metrics,
        "trained_at": datetime.now(timezone.utc).isoformat(),
        "software": {"name": "ViroFlow", "version": __version__, "author": __author__},
    }
    model_path.parent.mkdir(parents=True, exist_ok=True)
    sk["joblib"].dump(bundle, model_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps({key: value for key, value in bundle.items() if key != "model"}, indent=2) + "\n")
    return metrics


def predict_classifier(
    model_path: Path,
    features_path: Path,
    output_path: Path,
    sample_column: str = "sample",
) -> Path:
    sk = _require_sklearn()
    bundle = sk["joblib"].load(model_path)
    sample_ids, feature_columns, features, _ = read_feature_table(
        features_path,
        label_column=None,
        sample_column=sample_column,
        exclude_columns={bundle.get("label_column", "")},
    )
    expected = bundle["feature_columns"]
    if feature_columns != expected:
        raise MLError("feature schema mismatch; expected columns: " + ", ".join(expected))
    x_matrix = sk["np"].asarray(features, dtype=float)
    model = bundle["model"]
    probabilities = model.predict_proba(x_matrix)
    classes = [str(item) for item in model.classes_]
    predictions = model.predict(x_matrix)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["sample", "predicted_label", "confidence"] + [f"probability_{label}" for label in classes]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for sample, prediction, probability in zip(sample_ids, predictions, probabilities):
            row = {
                "sample": sample,
                "predicted_label": prediction,
                "confidence": round(float(max(probability)), 6),
            }
            row.update({f"probability_{label}": round(float(value), 6) for label, value in zip(classes, probability)})
            writer.writerow(row)
    return output_path


def score_anomalies(
    features_path: Path,
    output_path: Path,
    contamination: str | float = "auto",
    sample_column: str = "sample",
    exclude_columns: set[str] | None = None,
) -> Path:
    sk = _require_sklearn()
    metadata_columns = {"label", "labels", "class", "target", "outcome"}
    if exclude_columns:
        metadata_columns.update(exclude_columns)
    sample_ids, feature_columns, features, _ = read_feature_table(
        features_path,
        label_column=None,
        sample_column=sample_column,
        exclude_columns=metadata_columns,
    )
    x_matrix = sk["np"].asarray(features, dtype=float)
    model = sk["Pipeline"](
        [
            ("imputer", sk["SimpleImputer"]()),
            (
                "isolation_forest",
                sk["IsolationForest"](
                    n_estimators=500,
                    contamination=contamination,
                    random_state=42,
                ),
            ),
        ]
    )
    model.fit(x_matrix)
    decision = model.decision_function(x_matrix)
    raw_score = model.score_samples(x_matrix)
    flags = model.predict(x_matrix)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "anomaly_flag", "decision_score", "raw_score"],
            lineterminator="\n",
        )
        writer.writeheader()
        for sample, flag, decision_score, score in zip(sample_ids, flags, decision, raw_score):
            writer.writerow(
                {
                    "sample": sample,
                    "anomaly_flag": int(flag == -1),
                    "decision_score": round(float(decision_score), 6),
                    "raw_score": round(float(score), 6),
                }
            )
    metadata_path = output_path.with_suffix(output_path.suffix + ".metadata.json")
    metadata_path.write_text(
        json.dumps(
            {
                "method": "IsolationForest",
                "features": feature_columns,
                "contamination": contamination,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "software": {"name": "ViroFlow", "version": __version__, "author": __author__},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_path
