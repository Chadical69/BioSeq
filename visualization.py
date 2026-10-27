"""Plotly figures used by the BioSeq dashboard."""
from __future__ import annotations

from typing import Any

import plotly.graph_objects as go

PALETTE = {
    "A": "#14B8A6",
    "T": "#60A5FA",
    "G": "#A78BFA",
    "C": "#F59E0B",
    "GC": "#14B8A6",
    "AT": "#60A5FA",
}
LAYOUT = {
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
    "font": {"family": "Inter, ui-sans-serif, system-ui, sans-serif", "color": "#D9E4F2"},
    "margin": {"l": 24, "r": 24, "t": 48, "b": 28},
    "title": {"x": 0.02, "xanchor": "left", "font": {"size": 16}},
    "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
}


def _base_figure(fig: go.Figure) -> go.Figure:
    fig.update_layout(**LAYOUT)
    fig.update_xaxes(showgrid=False, zeroline=False, linecolor="rgba(148,163,184,.25)")
    fig.update_yaxes(gridcolor="rgba(148,163,184,.12)", zeroline=False, linecolor="rgba(148,163,184,.25)")
    return fig


def nucleotide_bar(counts: dict[str, int]) -> go.Figure:
    fig = go.Figure(go.Bar(
        x=list(counts.keys()), y=list(counts.values()),
        marker_color=[PALETTE[b] for b in counts],
        text=list(counts.values()), textposition="outside",
        hovertemplate="Base: %{x}<br>Count: %{y}<extra></extra>",
    ))
    fig.update_layout(title="Nucleotide counts", xaxis_title="Nucleotide", yaxis_title="Count", showlegend=False)
    return _base_figure(fig)


def nucleotide_donut(counts: dict[str, int]) -> go.Figure:
    labels = [base for base, value in counts.items() if value > 0]
    values = [counts[base] for base in labels]
    if not labels:
        labels, values = ["No canonical bases"], [1]
    fig = go.Figure(go.Pie(
        labels=labels, values=values, hole=0.65,
        marker_colors=[PALETTE.get(label, "#64748B") for label in labels],
        textinfo="label+percent", sort=False,
        hovertemplate="%{label}<br>Count: %{value}<br>Share: %{percent}<extra></extra>",
    ))
    fig.update_layout(title="Nucleotide composition", showlegend=True)
    return _base_figure(fig)


def gc_at_chart(gc_content: float, at_content: float) -> go.Figure:
    fig = go.Figure(go.Bar(
        x=["GC content", "AT content"], y=[gc_content, at_content],
        marker_color=[PALETTE["GC"], PALETTE["AT"]],
        text=[f"{gc_content:.2f}%", f"{at_content:.2f}%"], textposition="outside",
        hovertemplate="%{x}: %{y:.2f}%<extra></extra>",
    ))
    fig.update_layout(title="GC vs AT content", yaxis_title="Percentage (%)", yaxis_range=[0, 110], showlegend=False)
    return _base_figure(fig)


def codon_frequency_chart(codon_counts: dict[str, int], top_n: int = 15) -> go.Figure:
    items = list(codon_counts.items())[:top_n]
    fig = go.Figure(go.Bar(
        x=[item[0] for item in items], y=[item[1] for item in items],
        marker_color="#60A5FA", text=[item[1] for item in items], textposition="outside",
        hovertemplate="Codon: %{x}<br>Occurrences: %{y}<extra></extra>",
    ))
    fig.update_layout(title=f"Most frequent codons (top {top_n})", xaxis_title="Codon", yaxis_title="Occurrences", showlegend=False)
    return _base_figure(fig)
