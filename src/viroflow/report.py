from __future__ import annotations

import csv
import html
import json
from dataclasses import asdict
from pathlib import Path

from .models import SampleResult


def write_reports(result: SampleResult, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = asdict(result)
    json_path = output_dir / "analysis.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    mutation_path = output_dir / "mutations.tsv"
    with mutation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "level", "region", "mutation", "classification", "evidence"])
        for segment, values in result.segments.items():
            for mutation in values["nucleotide_changes"]:
                writer.writerow([result.sample, "nucleotide", segment, mutation, "", ""])
        for protein, values in result.proteins.items():
            antigenic = set(values.antigenic_site_changes)
            marker_map = {item["mutation"]: item["evidence"] for item in values.matched_markers}
            for mutation in values.amino_acid_changes:
                classes = []
                if mutation in antigenic:
                    classes.append("configured_antigenic_site")
                if mutation in marker_map:
                    classes.append("configured_escape_marker")
                writer.writerow(
                    [
                        result.sample,
                        "amino_acid",
                        protein,
                        mutation,
                        ",".join(classes),
                        marker_map.get(mutation, ""),
                    ]
                )

    html_path = output_dir / "report.html"
    html_path.write_text(_render_html(result), encoding="utf-8")
    return {"json": json_path, "tsv": mutation_path, "html": html_path}


def _render_html(result: SampleResult) -> str:
    def badge(label: str, value: object) -> str:
        return (
            '<div class="card"><div class="label">'
            + html.escape(label)
            + '</div><div class="value">'
            + html.escape(str(value))
            + "</div></div>"
        )

    segment_rows = "".join(
        "<tr>"
        f"<td>{html.escape(segment)}</td>"
        f"<td>{values['reference_length']}</td>"
        f"<td>{values['query_length']}</td>"
        f"<td>{100 * values['identity']:.2f}%</td>"
        f"<td>{100 * values['n_content']:.2f}%</td>"
        f"<td>{len(values['nucleotide_changes'])}</td>"
        "</tr>"
        for segment, values in result.segments.items()
    )
    protein_rows = "".join(
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f"<td>{len(values.amino_acid_changes)}</td>"
        f"<td>{html.escape(', '.join(values.antigenic_site_changes) or '—')}</td>"
        f"<td>{html.escape(', '.join(item['mutation'] for item in values.matched_markers) or '—')}</td>"
        "</tr>"
        for name, values in result.proteins.items()
    )
    warning_items = "".join(f"<li>{html.escape(item)}</li>" for item in result.warnings)
    shift_text = html.escape(result.shift_screen["interpretation"])
    escape_text = html.escape(result.vaccine_escape_screen["interpretation"])
    raw_json = html.escape(json.dumps(asdict(result), indent=2))
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ViroFlow report — {html.escape(result.sample)}</title>
<style>
:root{{--ink:#132238;--muted:#64748b;--line:#dbe3ed;--blue:#075985;--bg:#f6f9fc}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1100px;margin:auto;padding:32px 20px 64px}}h1{{margin:0 0 4px;font-size:30px}}
.subtitle{{color:var(--muted);margin-bottom:24px}}.grid{{display:grid;
grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}}
.card,section{{background:#fff;border:1px solid var(--line);border-radius:12px;
box-shadow:0 2px 8px #0f172a0d}}.card{{padding:16px}}.label{{color:var(--muted);
font-size:12px;text-transform:uppercase;letter-spacing:.05em}}.value{{font-size:24px;
font-weight:700;color:var(--blue)}}section{{padding:20px;margin-top:18px}}h2{{font-size:19px;
margin:0 0 12px}}table{{border-collapse:collapse;width:100%}}th,td{{padding:9px 10px;
border-bottom:1px solid var(--line);text-align:left}}th{{font-size:12px;color:var(--muted);
text-transform:uppercase}}.notice{{border-left:4px solid #f59e0b;padding-left:12px}}
details pre{{overflow:auto;background:#0f172a;color:#dbeafe;padding:16px;border-radius:8px}}
</style></head><body><main>
<h1>ViroFlow analysis</h1><div class="subtitle">Sample: {html.escape(result.sample)}</div>
<div class="grid">
{badge("Drift evidence index", result.drift["evidence_index"])}
{badge("Escape priority", result.vaccine_escape_screen["priority_band"])}
{badge("Escape score", result.vaccine_escape_screen["priority_score"])}
{badge("Reassortment signal", "candidate" if result.shift_screen["candidate_reassortment_signal"] else "not detected")}
</div>
<section><h2>Genome quality and nucleotide changes</h2><table><thead><tr>
<th>Segment</th><th>Reference nt</th><th>Query nt</th><th>Identity</th>
<th>N content</th><th>Changes</th></tr></thead><tbody>{segment_rows}</tbody></table></section>
<section><h2>Protein screening</h2><table><thead><tr><th>Protein</th><th>AA changes</th>
<th>Antigenic-site changes</th><th>Evidence markers</th></tr></thead>
<tbody>{protein_rows}</tbody></table></section>
<section><h2>Genotype</h2><p>{html.escape(result.genotype["composite"] or "unassigned")}</p>
<p>{shift_text}</p></section>
<section class="notice"><h2>Interpretation boundary</h2><p>{escape_text}</p>
<p>Sequence screens require laboratory, phylogenetic, epidemiologic, and expert review.</p></section>
<section><h2>Warnings</h2><ul>{warning_items or "<li>None</li>"}</ul></section>
<section><details><summary>Machine-readable result</summary><pre>{raw_json}</pre></details></section>
</main></body></html>
"""

