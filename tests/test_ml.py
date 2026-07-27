import csv
import json
from pathlib import Path

import pytest

from viroflow.ml import extract_features, predict_classifier, score_anomalies, train_classifier

pytest.importorskip("sklearn")


def _write_feature_table(path: Path) -> None:
    rows = [
        {"sample": "s1", "drift_index": 2.0, "escape_score": 4.0, "mean_depth": 120, "label": "low"},
        {"sample": "s2", "drift_index": 3.0, "escape_score": 5.0, "mean_depth": 110, "label": "low"},
        {"sample": "s3", "drift_index": 4.0, "escape_score": 6.0, "mean_depth": 100, "label": "low"},
        {"sample": "s4", "drift_index": 60.0, "escape_score": 70.0, "mean_depth": 95, "label": "high"},
        {"sample": "s5", "drift_index": 65.0, "escape_score": 80.0, "mean_depth": 90, "label": "high"},
        {"sample": "s6", "drift_index": 70.0, "escape_score": 85.0, "mean_depth": 85, "label": "high"},
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_train_and_predict_classifier(tmp_path: Path):
    features = tmp_path / "features.csv"
    model = tmp_path / "model.joblib"
    report = tmp_path / "report.json"
    predictions = tmp_path / "predictions.csv"
    _write_feature_table(features)

    metrics = train_classifier(features, "label", model, report)
    assert model.exists()
    assert metrics["samples"] == 6
    assert metrics["features"] == 3

    predict_classifier(model, features, predictions)
    output = predictions.read_text(encoding="utf-8")
    assert "predicted_label" in output
    assert "probability_high" in output


def test_anomaly_scores_are_written(tmp_path: Path):
    features = tmp_path / "features.csv"
    scores = tmp_path / "scores.csv"
    _write_feature_table(features)

    score_anomalies(features, scores, contamination="auto")
    assert scores.exists()
    assert scores.with_suffix(".csv.metadata.json").exists()
    assert "anomaly_flag" in scores.read_text(encoding="utf-8")


def test_extract_features_from_viroflow_results(tmp_path: Path):
    sample_report = tmp_path / "results" / "sample01" / "05_report"
    sample_report.mkdir(parents=True)
    (sample_report / "analysis.json").write_text(
        json.dumps(
            {
                "segments": {
                    "genome": {
                        "identity": 0.99,
                        "n_content": 0.01,
                        "callable_sites": 29000,
                        "reference_length": 29903,
                        "nucleotide_changes": ["A10G"],
                        "inserted_bases": 0,
                    }
                },
                "proteins": {
                    "S": {
                        "amino_acid_changes": ["D614G"],
                        "antigenic_site_changes": ["D614G"],
                        "matched_markers": [{"mutation": "D614G", "weight": 1.5}],
                    }
                },
                "drift": {"evidence_index": 12.0, "amino_acid_divergence": 0.001},
                "vaccine_escape_screen": {"priority_score": 8.0},
                "shift_screen": {"candidate_reassortment_signal": False},
                "genotype": {"segment_calls": {"genome": {"lineage": "A", "confident": True}}},
            }
        ),
        encoding="utf-8",
    )
    (sample_report / "qc_summary.json").write_text(
        json.dumps(
            {
                "assembly": {"contigs": 2, "total_bases": 29903, "n50": 20000, "gc_fraction": 0.38},
                "coverage": {
                    "mean_depth": 100,
                    "median_depth": 95,
                    "min_depth": 10,
                    "breadth_at_1x": 1,
                    "breadth_at_10x": 1,
                    "breadth_at_100x": 0.5,
                    "uniformity_fraction_0_2x_to_2x_mean": 0.9,
                },
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "cohort_features.csv"
    extract_features(tmp_path / "results", output)
    text = output.read_text(encoding="utf-8")
    assert "drift_index" in text
    assert "sample01" in text
