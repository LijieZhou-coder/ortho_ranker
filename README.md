# Ortho Ranker

**A high-precision homolog screening and candidate prioritization tool for plant comparative genomics, integrating reciprocal BLAST evidence with biologically interpretable expression ranking.**

Ortho Ranker is designed for a common and difficult real-world problem in plant molecular biology: after homology screening, multiple candidate genes often remain, especially in recently duplicated families. In these cases, sequence similarity alone is often not enough to decide which candidate should be prioritized for experimental validation.

Ortho Ranker addresses this by preserving a strict sequence-based evidence chain while adding an optional **Expression Ranking** layer that re-orders candidate homologs according to target-tissue abundance and background specificity.

This project is intended not only as a utility, but also as a transparent and teachable methodological framework for student training, candidate prioritization, and reproducible comparative analysis.

---

## Highlights

- **Forward BLASTP-based homolog retrieval**
- **Coverage-aware first-tier candidate selection**
- **Reciprocal BLASTP support classification**
- **Optional expression-based prioritization**
- **Wide-table expression matrix support**
- **White-list sample column parsing**
- **Transparent handling of missing expression IDs**
- **Human-readable summary outputs**
- **CLI workflow suitable for reproducible benchmarking**

---

## Why Ortho Ranker

In many plant systems, especially those involving:
- duplicated transcription factor families,
- specialized metabolism,
- anthocyanin or pigment biosynthesis,
- floral organ specification,
- stress-responsive paralogs,

the "best hit by sequence" is not always the best gene to validate first.

For example:
- one homolog may look strong by sequence but be broadly expressed across all tissues;
- another may be slightly less decisive in reciprocal support, but highly enriched in the target tissue of interest.

In practice, researchers often need a tool that can prioritize genes using **both**:
1. **sequence-level support**, and  
2. **biological expression context**.

Ortho Ranker was built for exactly this situation.

---

## Method overview

Ortho Ranker follows a layered evidence strategy:

1. **Forward BLASTP** against a target proteome
2. **Coverage-aware preprocessing** of top hits
3. **Breakpoint-based first-tier candidate selection**
4. **Reciprocal BLASTP** against a reference proteome
5. **Reciprocal support classification**
6. **Optional expression ranking** for candidate prioritization

The key principle is:

> **Expression does not redefine orthology.**  
> It acts only as a post hoc prioritization layer added after sequence-based candidate identification.

This preserves the original evidence chain while improving biological relevance in downstream ranking.

---

## Expression Ranking strategy

The current release supports a practical and interpretable ranking mode called **`combined`**.

Given:
- `target_mean` = mean expression across user-defined target samples
- `background_mean` = mean expression across user-defined background samples
- `pseudocount` = default 1.0

the score is:

\[
\text{score}
=
\log_2(\text{target\_mean} + p)
+
\log_2\left(
\frac{\text{target\_mean} + p}{\text{background\_mean} + p}
\right)
\]

This score rewards:
- **high abundance in the target tissue**
- **strong specificity relative to background tissues**

This is especially useful for avoiding two common mistakes:
- prioritizing broadly expressed housekeeping-like homologs;
- overvaluing extremely low-expression genes with misleadingly large fold changes.

### Current input metric

The current implementation supports **FPKM** as the expression input metric.

This should be interpreted as:
- a **candidate prioritization score**
- **not** a formal differential expression significance test

---

## Workflow overview

```text
Query protein(s)
    ↓
Forward BLASTP against target proteome
    ↓
Coverage-aware top-hit preprocessing
    ↓
First-tier candidate selection
    ↓
Export candidate FASTA
    ↓
Reciprocal BLASTP against reference proteome
    ↓
Reciprocal support summary
    ↓
(Optional) Expression Ranking
    ↓
Final candidate prioritization
```

---

## Installation

### Requirements

- Python 3.10+
- BLAST+ installed and available in `PATH`

### Clone the repository

```bash
git clone https://github.com/<your-username>/ortho_ranker.git
cd ortho_ranker
```

### Create and activate an environment

Using `conda`:

```bash
conda create -n ortho_ranker python=3.11 -y
conda activate ortho_ranker
```

Or using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate
```

### Install Ortho Ranker

#### Standard installation

Install the package and its runtime dependencies:

```bash
pip install -e .
```

This installs:
- `ortho_ranker`
- `pandas`
- `typer`
- `pyyaml`
- `biopython`

#### Development installation

Install the package together with development/test dependencies:

```bash
pip install -e .[dev]
```

This is recommended if you want to run the test suite or contribute to development.

### Install BLAST+

Ortho Ranker depends on BLAST+ for sequence similarity searches. If BLAST+ is not already installed, one convenient option is:

```bash
conda install -c bioconda blast -y
```

After installation, confirm that BLAST is available:

```bash
blastp -version
```

### Verify the installation

Check that the CLI entry point is available:

```bash
ortho-ranker --help
```

If you installed the development dependencies, you can also run:

```bash
pytest -q
```

---

## Quick Start

### 1. Prepare BLAST databases

Before running the pipeline, build BLAST databases for the target and reference proteomes if they do not already exist.

Example:

```bash
makeblastdb -in target_proteome.fa -dbtype prot -out test_data/blastdb/target_test_db
```

### 2. Validate configuration and paths

```bash
ortho-ranker assess config.test.yaml
```

This checks:
- config structure
- required files
- output directory
- BLAST database availability

### 3. Run the pipeline

```bash
ortho-ranker run config.test.yaml
```

If expression ranking is enabled, the workflow will automatically continue from reciprocal support into expression-based candidate ranking.

---

## Configuration example

A minimal expression section looks like this:

```yaml
expression:
  enabled: true
  matrix_file: test_data/expression/expression_test.tsv
  gene_id_column: gene_id
  expression_metric: FPKM
  target_samples:
    - petal_1
    - petal_2
  background_samples:
    - leaf_1
    - leaf_2
  ranking_strategy: combined
  min_expression_threshold: 0.0
  include_no_valid_hit: true
  pseudocount: 1.0
```

### Important notes

- `matrix_file` should be a **wide-format expression matrix**
- only configured sample columns are used
- extra columns such as `GO_annotation`, `KEGG`, or `COG_class` are ignored
- missing candidate IDs are retained and explicitly marked rather than silently removed

---

## Input expectations

### Sequence inputs

- query protein FASTA
- target proteome FASTA
- reference proteome FASTA
- corresponding BLAST databases

### Expression matrix

The current implementation expects a wide-format table such as:

```tsv
gene_id	petal_1	petal_2	leaf_1	leaf_2	GO_annotation	KEGG
target_hit_1	100	110	5	4	chlorophyll biosynthesis	ko00195
target_hit_2	60	55	20	18	stress response	ko04075
```

Only these columns are used for ranking:
- `gene_id_column`
- `target_samples`
- `background_samples`

Everything else is ignored.

---

## Main outputs

Typical output files include:

- `forward_blast.tsv`
- `forward_top_hits.tsv`
- `forward_candidates.tsv`
- `forward_breakpoint_summary.txt`
- `reciprocal_query_candidates.fa`
- `reciprocal_blast.tsv`
- `reciprocal_support.tsv`
- `reciprocal_summary.txt`

If expression ranking is enabled, two additional files are produced:

- `expression_ranking.tsv`
- `expression_summary.txt`

### `expression_ranking.tsv`

This table contains fields such as:
- `candidate_id`
- `expected_query_id`
- `support_level`
- `expression_data_found`
- `target_mean_fpkm`
- `background_mean_fpkm`
- `specificity_ratio`
- `expression_priority_score`
- `expression_rank`
- `recommended_by_expression`

### `expression_summary.txt`

This text summary reports:
- expression metric
- ranking strategy
- target/background sample names
- number of ranked candidates
- number of candidates with or without expression data
- top-ranked candidate
- top candidate target mean
- top candidate background mean
- top specificity ratio
- top expression score
- warnings when applicable

---

## Benchmark case

This repository includes a small teaching-oriented benchmark that illustrates why expression-aware prioritization is useful after reciprocal homolog screening.

### Biological scenario

After reciprocal analysis, two candidate homologs remain:

- **`target_hit_1`**: strongly enriched in the target tissue (petal), low in background tissues
- **`target_hit_2`**: expressed in the target tissue, but also broadly expressed in background tissues

Example summary:

```text
target_hit_1:
  target mean = 105.0
  background mean = 4.5
  specificity ratio ≈ 19.27

target_hit_2:
  target mean = 57.5
  background mean = 19.0
  specificity ratio ≈ 2.93
```

Expected biological interpretation:
- `target_hit_1` should be prioritized for validation
- `target_hit_2` remains a plausible homolog, but is less tissue-specific

### Pseudo-specific low-expression control

The test suite also includes a synthetic low-expression control:

- **`target_hit_3`**: very low abundance, superficially specific due to near-zero background

This is used to demonstrate an important principle:

> a near-zero-expression gene should not outrank a genuinely high-expression, target-enriched candidate simply because its background is also near zero.

This case helps justify the use of the `combined` score instead of naïve fold-change-only ranking.

---

## Testing

Run the full test suite with:

```bash
pytest -q
```

The current frozen benchmark/test suite covers:

- end-to-end CLI execution
- expression output generation
- explainability summary assertions
- missing-expression candidate retention
- clean failure on missing required expression columns
- low-expression pseudo-specific candidate behavior
- empty-candidate edge-case handling

This benchmark is intended to function as a compact, interpretable, teaching-grade mini-example rather than a large-scale statistical benchmark collection.

---

## Interpretation guide

### Reciprocal support labels

- **`strong`**: reciprocal top hit matches expectation and is sufficiently separated from the next hit
- **`ambiguous`**: reciprocal evidence exists, but not strongly enough to be decisive
- **`no_valid_hit`**: no reciprocal hit passes configured coverage thresholds

### Expression ranking output

Expression ranking should be interpreted as:
- a **priority recommendation**
- not an orthology redefinition
- not a formal DE significance result

This means a candidate with weaker expression context should not disappear; it should remain visible in the result table together with its support label and ranking evidence.

---

## Limitations

- The current release focuses on **protein-based homolog screening**
- Expression ranking currently supports **FPKM** input only
- The score is designed for **candidate prioritization**, not formal cross-sample statistical inference
- Expression ranking depends on correctly specified sample column names
- Reciprocal support and expression evidence should be interpreted together

---

## Citation and reuse

If you use Ortho Ranker in a manuscript, thesis, teaching material, or benchmark demonstration, please cite the repository URL and specify the version used.

Suggested citation format:

```text
Zhou L. Ortho Ranker: a high-precision homolog screening and candidate prioritization tool integrating reciprocal BLAST evidence with expression ranking. GitHub repository, 2026. <repo-url>
```

Please also cite the original tools and resources underlying your analysis where appropriate, including:
- BLAST+
- the transcriptome quantification method used to generate the expression matrix
- proteome and annotation sources used in the workflow

---

## Reproducibility notes

For reproducible use:
- keep config files under version control
- archive exact FASTA and expression input files
- record BLAST database versions
- preserve summary files together with ranking tables
- avoid manual editing of intermediate outputs

---

## Roadmap

Possible future directions include:

- TPM support
- normalized count support
- alternative ranking strategies
- richer report generation
- optional visualization modules
- larger benchmark collections
- tighter domain-architecture integration

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Acknowledgments

Ortho Ranker was shaped by practical needs in plant molecular biology, especially the recurring challenge of prioritizing duplicated homolog candidates for downstream functional validation.

The project is intended not only as a command-line tool, but also as a transparent and teachable methodological framework.

---

## Contact

Questions, suggestions, and benchmark contributions are welcome via GitHub Issues or Pull Requests.
