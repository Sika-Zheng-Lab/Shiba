"""A contrast-aware report for regression results (no pairwise PSI assumptions)."""
import html
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def write_report(input_dir, output_dir):
    import matplotlib.pyplot as plt
    import plotly.express as px

    source = Path(input_dir) / "splicing"
    output = Path(output_dir)
    (output / "png").mkdir(parents=True, exist_ok=True)
    (output / "pdf").mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source / "analysis.json").read_text())
    results = pd.read_csv(source / "event_statistics.tsv", sep="\t")
    summary = pd.read_csv(source / "summary.txt", sep="\t")
    predictions = pd.read_csv(source / "model_predictions.tsv", sep="\t")
    sections = []
    for i, (contrast, table) in enumerate(results.groupby("contrast", sort=False)):
        table = table.copy()
        table["-log10(q_event)"] = -np.log10(table.q_event.clip(lower=np.finfo(float).tiny))
        figure = px.scatter(table, x="conservative_component_effect", y="-log10(q_event)",
                            color="Diff events", symbol="event_type",
                            hover_data=["event_id", "estimate_min", "estimate_max", "q_event"],
                            labels={"conservative_component_effect": "Smallest absolute component effect (logit contrast)"},
                            title=str(contrast))
        sections.append(figure.to_html(full_html=False, include_plotlyjs=True if i == 0 else False))
    counts = results[results["Diff events"] == "Yes"].groupby(["contrast", "event_type"]).size()
    fig, ax = plt.subplots(figsize=(10, 5))
    if len(counts):
        counts.unstack(fill_value=0).plot.bar(ax=ax)
    else:
        ax.text(.5, .5, "No significant events", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("Significant events")
    ax.set_xlabel("Logit-scale contrast")
    fig.tight_layout()
    fig.savefig(output / "png" / "barplot_splicing_summary.png")
    fig.savefig(output / "pdf" / "barplot_splicing_summary.pdf")
    plt.close(fig)
    document = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>Shiba beta-binomial regression</title><style>
body{font:16px system-ui,sans-serif;margin:2rem;max-width:1400px}table{border-collapse:collapse}
th,td{padding:.4rem;border:1px solid #ddd}pre{white-space:pre-wrap}.table{overflow:auto}
</style><h1>Shiba beta-binomial regression</h1>"""
    document += "<p>Formula: <code>" + html.escape(manifest["formula"]) + "</code></p>"
    document += ("<p>P values test logit-scale contrasts. Model PSI summaries are equal-weight component means, "
                 "distinct from observed Shiba PSI. Missing/failed components prevent an event call. "
                 "FDR is adjusted across coverage-passing events from all event types separately for each contrast. "
                 "Events marked filtered_low_coverage are excluded from fitting and FDR correction.</p>")
    document += "<h2>Event summary</h2>" + summary.to_html(index=False, escape=True)
    document += "".join(sections)
    for title, filename, table in [("Event results", "event_statistics.tsv", results),
                                   ("Model PSI predictions and effects", "model_predictions.tsv", predictions)]:
        link = html.escape(os.path.relpath(source / filename, output), quote=True)
        document += f'<h2>{title}</h2><p><a href="{link}">Complete TSV table</a>. Showing up to 2,000 rows.</p>'
        document += "<div class='table'>" + table.head(2000).to_html(index=False, escape=True) + "</div>"
    document += "</html>"
    (output / "summary.html").write_text(document)
