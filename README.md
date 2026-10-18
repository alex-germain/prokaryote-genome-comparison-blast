# Prokaryotic genome comparison

A desktop application (Python + Tkinter) that compares two prokaryotic genomes
from NCBI RefSeq: it downloads them, aligns their proteomes with BLASTp,
displays the result as an interactive dotplot, annotates each protein with its
COG functional category, and exports a report as CSV or PDF.

Written as a university project (MSc Bioinformatics, Sorbonne University),
then cleaned up for publication.

![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## What it does

Given two species picked from a searchable list of every bacterial and archaeal
assembly in RefSeq (~400 000 entries), the application:

1. **Downloads** the protein set (`.faa`) and the gene annotations (`.gff`)
   from the NCBI FTP server.
2. **Runs BLASTp** of genome 1 against genome 2, with a user-adjustable
   e-value threshold.
3. **Draws a dotplot** — each point is a protein pair, placed at the real
   chromosomal coordinates of both proteins, so conserved gene order appears
   as a diagonal. Points are hoverable: the identity, e-value and COG category
   of a hit are shown on mouse-over.
4. **Detects synteny blocks** — runs of hits close to a common diagonal,
   highlighted in red.
5. **Assigns COG categories** by running RPS-BLAST against the NCBI COG
   database, and recolours the dotplot by functional category.
6. **Reports statistics** — mean, minimum and maximum percent identity, number
   of alignments, number of conserved proteins on each side, plus COG
   distribution pie charts for each genome and for the shared proteins.
7. **Exports** either a CSV table of every hit with its COG annotations, or an
   illustrated PDF report. The report includes an automatic reading of the COG
   pie charts: the most conserved functional categories among shared proteins,
   the category most unbalanced between the two proteomes, and the share of
   each proteome that actually received an annotation.

### Reading the output

Mean identity is turned into a plain-language statement: above 90% the two
genomes are most likely strains of the same species, 70–90% species of the same
genus, 50–70% different genera, below 50% different families. This is a rule of
thumb for orientation, not a taxonomic assignment.

The COG colour code is shared by the dotplot, the legend window, the pie charts
and the PDF, and the 26 colours are grouped by hue so that the broad functional
group of a point can be read without consulting the legend: blues for
information storage and processing, greens and teals for cellular processes and
signalling, oranges and reds for metabolism, purples for poorly characterised
categories.

---

## Requirements

**Python 3.8+** with:

```
requests
matplotlib
mplcursors
reportlab
```

`tkinter` is part of the standard library but is packaged separately on some
Linux distributions (`sudo apt install python3-tk` on Debian/Ubuntu).

**NCBI BLAST+** must be on your `PATH`, providing `makeblastdb`, `blastp` and
`rpsblast`:

```bash
# Debian / Ubuntu
sudo apt install ncbi-blast+

# conda
conda install -c bioconda blast

# Windows / macOS
# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

> On Debian and Ubuntu the package installs RPS-BLAST as `rpsblast+`, not
> `rpsblast`, because the name `rpsblast` belonged to the older `blast`
> package. The application looks for both names, so no manual symlink is
> needed.

**COG database** (only required for the functional analysis), from
<https://ftp.ncbi.nlm.nih.gov/pub/COG/>:

- `cognames2003-2014.tab` and `fun2003-2014.tab` in the project folder
- the RPS-BLAST COG profile database in `COG/`, named `Cog`

---

## Installation

```bash
git clone https://github.com/alex-germain/prokaryote-genome-comparison-blast.git
cd prokaryote-genome-comparison-blast

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python PROJET_PYTHON_GERMAIN.py
```

On the first launch the application downloads the RefSeq metadata for bacteria
and archaea (~200 MB, a few minutes) and builds a local index in `Local_DB/`.
Every later launch reads that local copy and starts immediately.

---

## Usage

1. **Select two genomes.** Optionally filter by domain, type a species name to
   narrow the list, then click an entry. A green background confirms the
   selection.
2. **Download the genomes.** Files land in `genomes/`.
3. **Run BLAST.** Adjust the e-value if needed (default `1e-20`). BLAST runs in
   a background thread, so the interface stays responsive.
4. **Visualise.** Simple dotplot, synteny regions, or COG colour code. The
   minimum-identity filter redraws the current plot with only the hits above
   the threshold.
5. **More info.** Statistics, automatic interpretation and the three COG pie
   charts.
6. **Save results** as CSV or PDF.

Results are cached: `genomes/`, `blast_results/` and `cog_results/` are reused
instead of being recomputed. Delete a file to force it to be regenerated.

---

## Project layout

| File | Role |
|---|---|
| `PROJET_PYTHON_GERMAIN.py` | Tkinter interface and application logic |
| `create_local_DB.py` | Builds the local RefSeq metadata index |
| `telechargement_genomes.py` | Downloads and decompresses `.faa` and `.gff` |
| `lancer_blast.py` | `makeblastdb` + `blastp` |
| `cog_analysis.py` | RPS-BLAST against the COG database, result parsing |
| `dotplot.py` | Dotplot coordinates and synteny block detection |
| `genome_stats.py` | BLAST statistics and COG distributions |
| `cog_style.py` | Single source of truth for COG colours, order and labels |
| `export_results.py` | CSV and PDF export |
| `tools.py` | Shared helpers: tool lookup, safe filenames, atomic downloads |

---

## Design notes

A few decisions are worth stating, because they are the ones that took the most
time to get right.

**Identity is averaged over best hits, not over all alignments.** A protein with
many paralogues would otherwise contribute as many values as it has hits and
dominate the mean. Each query protein contributes exactly one value, from its
highest-bitscore alignment.

**Three distinct counts are reported separately.** The number of alignments,
the number of conserved proteins in genome 1 and the number in genome 2 are
different quantities; conflating them makes a comparison look more symmetric
than it is.

**Every expensive result is cached, so every expensive result is written
atomically.** The application decides whether to redo work by testing whether
the output file exists. That test is only sound if a file never appears at its
final path until it is complete — otherwise an interrupted download leaves a
truncated file that all later runs treat as valid, and the error only surfaces
much later as an unreadable archive or a truncated FASTA. Downloads and
decompressions therefore write to a `.part` file and are renamed into place
with `os.replace()`, which is atomic. For the same reason an empty BLAST output
is treated as missing rather than as an empty result.

**External tool names are resolved, not hard-coded**, through
`tools.find_blast_tool()`, which tries every known name for a tool and raises an
error naming what to install rather than a bare `FileNotFoundError`. The
`rpsblast` / `rpsblast+` split above is the reason.

**Organism names are sanitised before being used in filenames.** NCBI names
routinely contain parentheses, apostrophes, colons and brackets; `:` and `?`
are legal on Linux and illegal on Windows, so a name that works on one machine
breaks on the other.

## Known limitations

- BLASTp of two complete proteomes takes several minutes on a laptop; this is
  BLAST's cost, not the application's.
- Synteny detection is a simple diagonal-binning heuristic, not a published
  synteny algorithm — it is meant to make conserved gene order visible, not to
  produce publication-grade synteny blocks.
- Only the single best COG hit per protein is kept, so multi-domain proteins
  are assigned one category.
- The automatic interpretations are descriptive summaries of the computed
  numbers, not statistical tests. A difference of a few percentage points
  between two COG profiles is reported as such, without a significance claim.
- The interface and reports are in English; inline code comments are in French.

---

## Author

Alexandre Germain — MSc Bioinformatics and Modelling, Sorbonne University.

Developed as a coursework project and subsequently refactored with AI
assistance (portability fixes, atomic downloads, unified COG colour handling).

Released under the MIT licence - see `LICENSE`.
