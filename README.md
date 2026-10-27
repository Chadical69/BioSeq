# BioSeq — DNA Sequence Analysis Toolkit

BioSeq is a local, interactive bioinformatics learning app for analyzing DNA sequences. It is built with Python, Streamlit, Biopython, Plotly, and Pandas.

## Features

- Paste plain DNA or FASTA text, or upload `.fasta`, `.fa`, `.fna`, and `.txt` files
- Validate standard and IUPAC ambiguous DNA symbols
- Calculate sequence length, nucleotide counts, nucleotide percentages, GC content, AT content, and ambiguous-base count
- Generate complementary DNA and reverse-complement sequences
- Transcribe DNA to RNA using a coding-strand or template-strand convention
- Translate a DNA coding sequence in reading frames +1, +2, or +3
- Find ATG start-codon positions and standard stop-codon positions in the selected frame
- Show an interactive codon frequency table
- View Plotly charts for nucleotide counts, composition, GC/AT content, and codon frequencies
- Download CSV, TXT, FASTA, and JSON files
- Load included sample sequences for demonstration

## Requirements

- Windows 10/11, macOS, or Linux
- Python 3.10 or newer recommended
- Internet is not required after Python and the dependencies have been installed. The UI imports Google Fonts for typography; if offline, it falls back to system fonts.

## Windows quick start

1. Extract the `BioSeq` folder to a location you can write to.
2. Install Python from <https://www.python.org/downloads/> if it is not already installed. During setup, enable **Add Python to PATH**.
3. Double-click `install_and_run.bat`.
4. The script creates `.venv`, installs dependencies, and starts the app in your default browser.

You can also run it manually from PowerShell or Command Prompt:

```powershell
cd BioSeq
py -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -m streamlit run app.py
```

If the `py` command is unavailable, use `python` instead.

## macOS / Linux quick start

```bash
cd BioSeq
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -m streamlit run app.py
```

## Input notes

- Plain sequences and FASTA-formatted sequences are accepted.
- FASTA headers (lines beginning with `>`) and FASTA comment lines (beginning with `;`) are ignored during normalization.
- Whitespace is removed and letters are uppercased. Invalid characters are reported rather than silently discarded.
- Supported symbols are `A C G T R Y S W K M B D H V N`.
- GC and AT percentages use only unambiguous A/T/G/C bases as the denominator; ambiguous symbols are counted separately.
- Translation assumes the entered sequence is a coding strand in the 5′→3′ direction. The transcription selector independently controls how RNA is produced.

## Project structure

```text
BioSeq/
├── app.py
├── dna_analysis.py
├── visualization.py
├── export_utils.py
├── requirements.txt
├── install_and_run.bat
├── README.md
├── sample_sequences/
│   └── sample_dna.fasta
└── tests/
    └── test_dna_analysis.py
```

## Run the tests

After installing requirements, run:

```bash
python -m unittest discover -s tests -v
```

The tests cover input normalization, validation, base statistics, IUPAC complements, transcription, and translation.

## Scientific limitations

BioSeq is an educational sequence utility, not a clinical or diagnostic tool. It does not predict gene expression or establish that a translated sequence is a functional protein. A start codon or stop codon occurring in a chosen frame is not sufficient to identify a complete gene.
