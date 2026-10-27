from __future__ import annotations

from pathlib import Path
import json

import pandas as pd
import streamlit as st

from dna_analysis import analyze_sequence, normalize_sequence, validate_sequence
from export_utils import build_csv_report, build_fasta, build_json_report, build_text_report
from visualization import codon_frequency_chart, gc_at_chart, nucleotide_bar, nucleotide_donut

APP_DIR = Path(__file__).resolve().parent
DEFAULT_DNA = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
SAMPLE_SEQUENCES = {
    "Short coding sequence": "ATGGCTTAA",
    "Demo sequence (recommended)": DEFAULT_DNA,
    "Longer mixed sequence": "ATGCGTACCGATGCTAGCTAGGCTAATCGGATCCGATGCGTAA",
    "Ambiguous bases example": "ATGGCNTTRGCTAA",
}
STRAND_OPTIONS = ["Coding strand (5′→3′)", "Template strand (3′→5′)"]

st.set_page_config(
    page_title="BioSeq | DNA Sequence Workbench",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Manrope:wght@400;500;600;700;800&display=swap');
    :root { --bio-ink:#eaf2ff; --bio-muted:#9aabc5; --bio-teal:#2dd4bf; --bio-blue:#60a5fa; }
    .stApp { background: radial-gradient(ellipse at 14% 0%, rgba(20,184,166,.12), transparent 36%), #08111f; color:var(--bio-ink); }
    [data-testid="stHeader"] { background:rgba(8,17,31,.80); }
    [data-testid="stSidebar"] { background:linear-gradient(180deg,#0d1a2b,#0a1422); border-right:1px solid rgba(148,163,184,.12); }
    .block-container { padding-top:1.5rem; padding-bottom:3rem; max-width:1440px; }
    h1,h2,h3,h4,p,li,label { font-family:'Manrope',sans-serif; }
    .hero { position:relative; overflow:hidden; padding:1.6rem 1.8rem; border:1px solid rgba(45,212,191,.23); border-radius:22px; background:linear-gradient(120deg,rgba(15,39,58,.98),rgba(13,28,48,.95) 60%,rgba(15,56,62,.8)); margin-bottom:1.2rem; }
    .hero:after { content:'DNA'; position:absolute; right:26px; top:0px; font-size:5rem; font-weight:800; opacity:.055; letter-spacing:.12em; color:#fff; pointer-events:none; }
    .hero .eyebrow { color:#5eead4; letter-spacing:.18em; font-size:.72rem; font-weight:800; margin-bottom:.45rem; }
    .hero h1 { color:#f4f8ff; margin:0; font-size:2.4rem; letter-spacing:-.06em; font-weight:800; }
    .hero p { color:#aabbd1; font-size:1rem; margin:.45rem 0 0; max-width:780px; line-height:1.65; }
    .section-kicker { color:#5eead4; text-transform:uppercase; letter-spacing:.14em; font-size:.7rem; font-weight:800; margin-bottom:.25rem; }
    .section-title { color:#eef5ff; font-size:1.25rem; font-weight:800; letter-spacing:-.025em; margin-bottom:.2rem; }
    .muted { color:#9aabc5; font-size:.9rem; line-height:1.55; }
    .metric-card { border:1px solid rgba(148,163,184,.16); border-radius:16px; background:linear-gradient(145deg,rgba(20,35,55,.95),rgba(12,25,41,.95)); padding:1rem 1.1rem; min-height:116px; }
    .metric-label { color:#9aabc5; font-size:.78rem; font-weight:700; margin-bottom:.55rem; }
    .metric-value { color:#f2f7ff; font-family:'Manrope',sans-serif; font-weight:800; font-size:1.85rem; letter-spacing:-.05em; line-height:1.15; }
    .metric-caption { color:#8ca0ba; font-size:.74rem; margin-top:.45rem; }
    .info-card { border:1px solid rgba(148,163,184,.15); border-radius:16px; padding:1rem 1.15rem; background:rgba(15,29,47,.72); margin:.25rem 0 .75rem; }
    .pill { display:inline-block; border-radius:999px; padding:.28rem .65rem; border:1px solid rgba(45,212,191,.25); background:rgba(20,184,166,.09); color:#73f1dc; font-size:.72rem; font-weight:700; margin-right:.35rem; }
    div[data-testid="stTabs"] button { font-weight:700; }
    div[data-testid="stTabs"] button[aria-selected="true"] { color:#5eead4; }
    div.stButton > button[kind="primary"], div.stFormSubmitButton > button { background:linear-gradient(100deg,#0d9488,#14b8a6); color:#031a1c; border:0; border-radius:11px; font-weight:800; min-height:2.8rem; }
    div.stButton > button, div[data-testid="stDownloadButton"] button { border-radius:10px; font-weight:700; }
    div[data-testid="stMetric"] { background:rgba(15,29,47,.8); padding:1rem; border:1px solid rgba(148,163,184,.16); border-radius:15px; }
    [data-testid="stDataFrame"] { border:1px solid rgba(148,163,184,.15); border-radius:12px; overflow:hidden; }
    .sequence-box { font-family:'DM Mono',monospace; font-size:.9rem; line-height:1.8; overflow-wrap:anywhere; padding:1rem; border:1px solid rgba(96,165,250,.2); background:rgba(5,13,25,.8); border-radius:12px; color:#c8e7ff; }
    .small-note { color:#91a4bf; font-size:.78rem; line-height:1.5; }
    hr { border-color:rgba(148,163,184,.13); }
    </style>
    """,
    unsafe_allow_html=True,
)

if "dna_input" not in st.session_state:
    st.session_state.dna_input = DEFAULT_DNA
if "active_raw" not in st.session_state:
    st.session_state.active_raw = DEFAULT_DNA
if "sample_choice" not in st.session_state:
    st.session_state.sample_choice = "Demo sequence (recommended)"


def _load_sample() -> None:
    selected = st.session_state.get("sample_choice", "Demo sequence (recommended)")
    sample = SAMPLE_SEQUENCES[selected]
    st.session_state.dna_input = sample
    st.session_state.active_raw = sample
    st.session_state.sample_loaded = selected

with st.sidebar:
    st.markdown("<div style='font-size:1.35rem;font-weight:800;color:#eaf2ff;'>🧬 BioSeq</div>", unsafe_allow_html=True)
    st.caption("DNA sequence analysis toolkit")
    st.markdown("---")
    st.markdown("**Quick start**")
    st.selectbox("Sample sequence", list(SAMPLE_SEQUENCES.keys()), key="sample_choice")
    st.button("Load selected sample", on_click=_load_sample, use_container_width=True)
    st.markdown("---")
    st.markdown("**Analysis settings**")
    strand_convention = st.selectbox(
        "DNA strand for transcription",
        STRAND_OPTIONS,
        help="Coding strand: enter the coding DNA strand 5′→3′. Template strand: enter the template strand 3′→5′.",
    )
    frame_label = st.selectbox(
        "Translation reading frame",
        ["+1 (start at base 1)", "+2 (start at base 2)", "+3 (start at base 3)"],
        help="Translation is calculated on the entered sequence as a coding strand in the 5′→3′ direction.",
    )
    reading_frame = int(frame_label[1])
    st.markdown("---")
    st.markdown("<div class='small-note'>Runs locally on your computer. No account or paid API required.</div>", unsafe_allow_html=True)

st.markdown(
    """
    <div class="hero">
      <div class="eyebrow">BIOINFORMATICS · SEQUENCE WORKBENCH</div>
      <h1>BioSeq</h1>
      <p>Turn a DNA string into meaningful biological summaries. Explore base composition, complementary strands, transcription, translation, interactive charts, and downloadable reports in one workspace.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.form("sequence_input_form", clear_on_submit=False):
    st.markdown("<div class='section-kicker'>Input</div><div class='section-title'>Bring a DNA sequence</div>", unsafe_allow_html=True)
    st.caption("Paste a raw DNA sequence or FASTA text, or upload a .fasta, .fa, or .txt file. IUPAC ambiguous symbols are supported.")
    typed_sequence = st.text_area(
        "DNA / FASTA sequence",
        key="dna_input",
        height=130,
        label_visibility="collapsed",
        placeholder=">my_sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
    )
    uploaded_file = st.file_uploader("Or upload a sequence file", type=["fasta", "fa", "fna", "txt"])
    analyze_clicked = st.form_submit_button("Analyze sequence", type="primary", use_container_width=True)

if analyze_clicked:
    if uploaded_file is not None:
        try:
            st.session_state.active_raw = uploaded_file.getvalue().decode("utf-8-sig")
        except UnicodeDecodeError:
            st.session_state.active_raw = uploaded_file.getvalue().decode("latin-1")
    else:
        st.session_state.active_raw = typed_sequence

raw_active = st.session_state.active_raw
sequence = normalize_sequence(raw_active)
validation_errors = validate_sequence(sequence)

if validation_errors:
    for err in validation_errors:
        st.error(err)
    st.info("Tip: use the sample selector in the sidebar to load a working example.")
    st.stop()

try:
    result = analyze_sequence(raw_active, strand_convention, reading_frame)
except ValueError as exc:
    st.error(str(exc))
    st.stop()

stats = result["statistics"]
translation = result["translation"]
seq = result["sequence"]

st.markdown(
    f"<span class='pill'>Sequence loaded</span><span class='pill'>{stats['sequence_length']:,} bases</span><span class='pill'>Frame +{reading_frame}</span>",
    unsafe_allow_html=True,
)
st.write("")

tabs = st.tabs(["Overview", "DNA operations", "Translation", "Visual analytics", "Export & help"])

with tabs[0]:
    st.markdown("<div class='section-kicker'>At a glance</div><div class='section-title'>Sequence overview</div>", unsafe_allow_html=True)
    st.markdown("<p class='muted'>A high-level snapshot of your current input. Percentages below use unambiguous A, T, G and C bases; ambiguous IUPAC bases are counted separately.</p>", unsafe_allow_html=True)
    metric_cols = st.columns(4)
    metric_data = [
        ("Sequence length", f"{stats['sequence_length']:,}", "Total nucleotide symbols"),
        ("GC content", f"{stats['gc_content']:.2f}%", "G + C among A/T/G/C bases"),
        ("AT content", f"{stats['at_content']:.2f}%", "A + T among A/T/G/C bases"),
        ("Ambiguous bases", f"{stats['ambiguous_bases']:,}", "IUPAC symbols other than A/T/G/C"),
    ]
    for col, (label, value, caption) in zip(metric_cols, metric_data):
        with col:
            st.markdown(
                f"<div class='metric-card'><div class='metric-label'>{label}</div><div class='metric-value'>{value}</div><div class='metric-caption'>{caption}</div></div>",
                unsafe_allow_html=True,
            )
    st.write("")
    left, right = st.columns([1.05, 1])
    with left:
        st.markdown("<div class='section-title'>Nucleotide counts</div>", unsafe_allow_html=True)
        st.dataframe(
            pd.DataFrame({"Base": list(stats["counts"].keys()), "Count": list(stats["counts"].values()), "Percent of canonical bases": [f"{stats['percentages'][b]:.2f}%" for b in stats["counts"]]}),
            hide_index=True, use_container_width=True,
        )
    with right:
        st.markdown("<div class='section-title'>Sequence preview</div>", unsafe_allow_html=True)
        preview = "\n".join(seq[i:i+60] for i in range(0, min(len(seq), 240), 60))
        if len(seq) > 240:
            preview += "\n… (preview limited to first 240 bases)"
        st.markdown(f"<div class='sequence-box'>{preview}</div>", unsafe_allow_html=True)
    st.markdown("<div class='info-card'><b>How to read this:</b> GC content describes the proportion of G and C among canonical bases. It is a descriptive sequence statistic, not by itself a diagnosis or proof of gene function.</div>", unsafe_allow_html=True)

with tabs[1]:
    st.markdown("<div class='section-kicker'>Sequence transformations</div><div class='section-title'>DNA, complement, and RNA</div>", unsafe_allow_html=True)
    st.caption("The complement is shown in the same display orientation as the input. Reverse complement reverses that complementary string.")
    op_tabs = st.tabs(["Original DNA", "Complement", "Reverse complement", "RNA transcript"])
    outputs = [
        (seq, "Input sequence in the entered order"),
        (result["complement"], "Complement using IUPAC base-pairing rules"),
        (result["reverse_complement"], "Reverse complement, often useful for sequence searches"),
        (result["rna_sequence"], f"RNA from {strand_convention}"),
    ]
    output_labels = ["Original DNA", "Complement", "Reverse complement", "RNA transcript"]
    for pane, (value, description), output_label in zip(op_tabs, outputs, output_labels):
        with pane:
            st.caption(description)
            st.code(value, language=None)
            slug = output_label.lower().replace(" ", "_")
            st.download_button(
                f"Download {output_label.lower()}",
                data=value + "\n",
                file_name=f"bioseq_{slug}.txt",
                mime="text/plain",
                key=f"download_seq_{slug}",
            )

with tabs[2]:
    st.markdown("<div class='section-kicker'>Genetic code</div><div class='section-title'>Translation results</div>", unsafe_allow_html=True)
    st.info("Translation treats the entered sequence as a DNA coding strand in the 5′→3′ direction. The transcription strand setting does not change translation.")
    trans_cols = st.columns(3)
    with trans_cols[0]:
        st.metric("Complete codons", f"{translation['codon_count']:,}")
    with trans_cols[1]:
        st.metric("Start codons (ATG)", str(len(translation["start_codon_positions"])))
    with trans_cols[2]:
        st.metric("Stop codons", str(len(translation["stop_codon_positions"])))
    st.markdown("**Amino acid sequence**")
    if translation["protein_sequence"]:
        st.code(translation["protein_sequence"], language=None)
    else:
        st.warning("This reading frame does not contain a complete codon.")
    if translation["trailing_bases_ignored"]:
        st.caption(f"{translation['trailing_bases_ignored']} trailing base(s) could not form a complete codon and were ignored.")
    start_text = ", ".join(map(str, translation["start_codon_positions"])) or "None found"
    stop_text = ", ".join(map(str, translation["stop_codon_positions"])) or "None found"
    left, right = st.columns(2)
    with left:
        st.markdown("**Start codon positions (1-based)**")
        st.write(start_text)
    with right:
        st.markdown("**Stop codon positions (1-based)**")
        st.write(stop_text)
    st.markdown("**Codon frequency table**")
    codon_df = pd.DataFrame([{"Codon": codon, "Occurrences": count} for codon, count in translation["codon_counts"].items()])
    if not codon_df.empty:
        st.dataframe(codon_df, hide_index=True, use_container_width=True)
    else:
        st.caption("No complete codons in this reading frame.")
    st.markdown("<div class='info-card'><b>Biological note:</b> A codon is a group of three nucleotides. The '*' symbol indicates a stop codon. Finding an ATG or stop codon does not by itself establish that the sequence is a complete, expressed gene.</div>", unsafe_allow_html=True)

with tabs[3]:
    st.markdown("<div class='section-kicker'>Interactive plots</div><div class='section-title'>Visual analytics</div>", unsafe_allow_html=True)
    st.caption("Charts update when you analyze a different sequence or change the reading frame.")
    chart_left, chart_right = st.columns(2)
    with chart_left:
        st.plotly_chart(nucleotide_bar(stats["counts"]), use_container_width=True)
    with chart_right:
        st.plotly_chart(nucleotide_donut(stats["counts"]), use_container_width=True)
    chart_left, chart_right = st.columns(2)
    with chart_left:
        st.plotly_chart(gc_at_chart(stats["gc_content"], stats["at_content"]), use_container_width=True)
    with chart_right:
        st.plotly_chart(codon_frequency_chart(translation["codon_counts"]), use_container_width=True)
    st.markdown("<div class='small-note'>Composition charts exclude ambiguous bases from the A/T/G/C percentage denominator. Codon frequencies are calculated from complete triplets in the selected reading frame.</div>", unsafe_allow_html=True)
with tabs[4]:
    st.markdown("<div class='section-kicker'>Take your results with you</div><div class='section-title'>Export analysis</div>", unsafe_allow_html=True)
    csv_report = build_csv_report(result)
    text_report = build_text_report(result)
    fasta_report = build_fasta(seq)
    json_report = build_json_report(result)
    exp_cols = st.columns(2)
    with exp_cols[0]:
        st.markdown("<div class='info-card'><b>CSV analysis report</b><br><span class='muted'>Statistics, nucleotide counts and percentages, sequence outputs, translation summary, and codon frequencies.</span></div>", unsafe_allow_html=True)
        st.download_button("⬇ Download CSV report", data=csv_report, file_name="bioseq_analysis.csv", mime="text/csv", use_container_width=True)
        st.markdown("<div class='info-card'><b>Plain-text report</b><br><span class='muted'>A readable summary suitable for a lab notebook or project appendix.</span></div>", unsafe_allow_html=True)
        st.download_button("⬇ Download TXT report", data=text_report, file_name="bioseq_report.txt", mime="text/plain", use_container_width=True)
    with exp_cols[1]:
        st.markdown("<div class='info-card'><b>FASTA sequence</b><br><span class='muted'>The cleaned DNA sequence with a FASTA header and wrapped sequence lines.</span></div>", unsafe_allow_html=True)
        st.download_button("⬇ Download FASTA", data=fasta_report, file_name="bioseq_sequence.fasta", mime="text/plain", use_container_width=True)
        st.markdown("<div class='info-card'><b>JSON analysis</b><br><span class='muted'>Structured output useful for reuse in other scripts or applications.</span></div>", unsafe_allow_html=True)
        st.download_button("⬇ Download JSON", data=json_report, file_name="bioseq_analysis.json", mime="application/json", use_container_width=True)
    st.markdown("---")
    with st.expander("Biology glossary"):
        st.markdown("""
        - **Nucleotide:** A DNA building block represented here by A, T, G, or C.
        - **GC content:** Percentage of unambiguous bases that are G or C.
        - **Complementary strand:** The sequence paired to the DNA input by A–T and G–C rules.
        - **Reverse complement:** The complementary strand reversed into the opposite orientation.
        - **Transcription:** Producing RNA from DNA; the convention depends on whether coding or template DNA is entered.
        - **Codon:** A three-base unit read during translation.
        - **Reading frame:** The position from which triplets are grouped (+1, +2, or +3).
        """)
    with st.expander("Assumptions and limitations"):
        st.markdown("""
        - Input is treated as a single DNA sequence, not an alignment or a full genome annotation.
        - IUPAC ambiguous symbols are accepted and preserved in sequence transformations. For composition percentages, only A/T/G/C count toward the denominator.
        - Translation uses the standard genetic code and the selected positive reading frame. It does not automatically identify the biologically correct open reading frame.
        - The default start-codon count reports ATG triplets in the selected frame; stop codons are TAA, TAG, and TGA.
        - Results are educational sequence analysis and are not medical or diagnostic interpretations.
        """)
    st.markdown("**Software:**Made with Python · Streamlit · Biopython · Plotly · Pandas")
st.markdown("---")
st.markdown("<div class='small-note' style='text-align:center'>BioSeq · DNA Sequence Analysis Toolkit ·</div>", unsafe_allow_html=True)
