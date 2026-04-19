from pathlib import Path
import subprocess
import pandas as pd


BLAST_OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def run_forward_blast(config: dict) -> dict:
    query_fasta = Path(config["paths"]["query_fasta"]).expanduser()
    target_db = Path(config["paths"]["target_blastdb_prefix"]).expanduser()
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    evalue_cutoff = config["thresholds"]["evalue_cutoff"]

    output_file = output_dir / "forward_blast.tsv"

    outfmt = "6 " + " ".join(BLAST_OUTFMT_FIELDS)

    cmd = [
        "blastp",
        "-query", str(query_fasta),
        "-db", str(target_db),
        "-evalue", str(evalue_cutoff),
        "-outfmt", outfmt,
        "-out", str(output_file),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    return {
        "command": " ".join(cmd),
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_file": str(output_file),
    }


def preprocess_forward_hits(config: dict) -> dict:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_file = output_dir / "forward_blast.tsv"
    top_hits_file = output_dir / "forward_top_hits.tsv"

    max_hits = int(config["thresholds"]["max_hits_to_keep"])

    if not raw_file.exists():
        raise FileNotFoundError(f"Forward BLAST output not found: {raw_file}")

    if raw_file.stat().st_size == 0:
        df = pd.DataFrame(columns=BLAST_OUTFMT_FIELDS + ["qcov", "scov", "rank"])
        df.to_csv(top_hits_file, sep="\t", index=False)
        return {
            "raw_rows": 0,
            "top_rows": 0,
            "top_hits_file": str(top_hits_file),
        }

    df = pd.read_csv(raw_file, sep="\t", header=None, names=BLAST_OUTFMT_FIELDS)

    numeric_cols = [
        "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "qlen", "slen"
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["qcov"] = df["length"] / df["qlen"]
    df["scov"] = df["length"] / df["slen"]

    df = df.sort_values(
        by=["qseqid", "bitscore", "evalue", "pident", "qcov", "scov"],
        ascending=[True, False, True, False, False, False]
    ).copy()

    df["rank"] = df.groupby("qseqid").cumcount() + 1
    top_df = df[df["rank"] <= max_hits].copy()

    top_df.to_csv(top_hits_file, sep="\t", index=False)

    return {
        "raw_rows": len(df),
        "top_rows": len(top_df),
        "top_hits_file": str(top_hits_file),
    }
