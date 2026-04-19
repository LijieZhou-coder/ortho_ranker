from pathlib import Path
import subprocess
import pandas as pd
from Bio import SeqIO


RECIPROCAL_OUTFMT_FIELDS = [
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


def export_candidate_fasta(config: dict) -> dict:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    target_proteome = Path(config["paths"]["target_proteome"]).expanduser()

    candidates_file = output_dir / "forward_candidates.tsv"
    reciprocal_fasta = output_dir / "reciprocal_query_candidates.fa"

    if not candidates_file.exists():
        raise FileNotFoundError(f"Candidate file not found: {candidates_file}")

    cand_df = pd.read_csv(candidates_file, sep="\t")
    selected_df = cand_df.loc[cand_df["selected_as_candidate"] == True].copy()
    selected_ids = set(selected_df["sseqid"].astype(str))

    records = []
    for record in SeqIO.parse(str(target_proteome), "fasta"):
        if record.id in selected_ids:
            records.append(record)

    SeqIO.write(records, str(reciprocal_fasta), "fasta")

    return {
        "candidate_count": len(selected_ids),
        "exported_count": len(records),
        "reciprocal_fasta": str(reciprocal_fasta),
    }


def run_reciprocal_blast(config: dict) -> dict:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    reciprocal_fasta = output_dir / "reciprocal_query_candidates.fa"
    reference_db = Path(config["paths"]["reference_blastdb_prefix"]).expanduser()

    evalue_cutoff = config["thresholds"]["evalue_cutoff"]
    output_file = output_dir / "reciprocal_blast.tsv"

    outfmt = "6 " + " ".join(RECIPROCAL_OUTFMT_FIELDS)

    cmd = [
        "blastp",
        "-query", str(reciprocal_fasta),
        "-db", str(reference_db),
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


def summarize_reciprocal_support(config: dict) -> dict:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    reciprocal_file = output_dir / "reciprocal_blast.tsv"
    support_file = output_dir / "reciprocal_support.tsv"
    summary_file = output_dir / "reciprocal_summary.txt"

    gap_cutoff = float(config["thresholds"]["reciprocal_top_gap_cutoff"])
    rqcov_cutoff = float(config["thresholds"]["reciprocal_query_cov_cutoff"])
    rscov_cutoff = float(config["thresholds"]["reciprocal_subject_cov_cutoff"])

    if not reciprocal_file.exists():
        raise FileNotFoundError(f"Reciprocal BLAST output not found: {reciprocal_file}")

    if reciprocal_file.stat().st_size == 0:
        empty_df = pd.DataFrame(columns=[
            "candidate_id", "expected_query_id", "top1_ref", "top2_ref",
            "top1_bitscore", "top2_bitscore", "gap_ratio",
            "top1_matches_expected", "support_level"
        ])
        empty_df.to_csv(support_file, sep="\t", index=False)
        with open(summary_file, "w", encoding="utf-8") as f:
            f.write("No reciprocal BLAST hits found.\n")
        return {
            "candidates_evaluated": 0,
            "strong_count": 0,
            "ambiguous_count": 0,
            "no_valid_hit_count": 0,
            "support_file": str(support_file),
            "summary_file": str(summary_file),
        }

    df = pd.read_csv(reciprocal_file, sep="\t", header=None, names=RECIPROCAL_OUTFMT_FIELDS)
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["length"] = pd.to_numeric(df["length"], errors="coerce")
    df["qlen"] = pd.to_numeric(df["qlen"], errors="coerce")
    df["slen"] = pd.to_numeric(df["slen"], errors="coerce")

    df["qcov"] = df["length"] / df["qlen"]
    df["scov"] = df["length"] / df["slen"]
    df["passes_reciprocal_cov_filter"] = (
        (df["qcov"] >= rqcov_cutoff) & (df["scov"] >= rscov_cutoff)
    )
    candidates_file = output_dir / "forward_candidates.tsv"
    cand_df = pd.read_csv(candidates_file, sep="\t")

    provenance_df = (
        cand_df.loc[cand_df["selected_as_candidate"] == True, ["sseqid", "expected_query_id"]]
        .drop_duplicates()
        .rename(columns={"sseqid": "candidate_id"})
    )

    df = df.sort_values(
        by=["qseqid", "bitscore", "evalue"],
        ascending=[True, False, True]
    ).copy()

    rows = []
    summary_lines = []

    for candidate_id, sub in df.groupby("qseqid", sort=False):
        sub = sub.reset_index(drop=True)

        match = provenance_df.loc[provenance_df["candidate_id"] == candidate_id, "expected_query_id"]

        if len(match) == 0:
            expected_query_id = pd.NA
        else:
            expected_query_id = match.iloc[0]

        filtered = sub[sub["passes_reciprocal_cov_filter"]].copy()
        filtered = filtered.sort_values(
            by=["bitscore", "evalue", "qcov", "scov"],
            ascending=[False, True, False, False]
        ).reset_index(drop=True)

        if filtered.empty:
            top1_ref = pd.NA
            top2_ref = pd.NA
            top1_bitscore = pd.NA
            top2_bitscore = pd.NA
            gap_ratio = pd.NA
            top1_matches_expected = False
            support_level = "no_valid_hit"
        else:
            top1_ref = filtered.loc[0, "sseqid"] if len(filtered) >= 1 else pd.NA
            top1_bitscore = filtered.loc[0, "bitscore"] if len(filtered) >= 1 else pd.NA

            top2_ref = filtered.loc[1, "sseqid"] if len(filtered) >= 2 else pd.NA
            top2_bitscore = filtered.loc[1, "bitscore"] if len(filtered) >= 2 else pd.NA

            if len(filtered) >= 2 and pd.notna(top1_bitscore) and top1_bitscore != 0:
                gap_ratio = (top1_bitscore - top2_bitscore) / top1_bitscore
            else:
                gap_ratio = pd.NA

            top1_matches_expected = (top1_ref == expected_query_id)

            if top1_matches_expected and pd.notna(gap_ratio) and gap_ratio >= gap_cutoff:
                support_level = "strong"
            else:
                support_level = "ambiguous"

        rows.append({
            "candidate_id": candidate_id,
            "expected_query_id": expected_query_id,
            "top1_ref": top1_ref,
            "top2_ref": top2_ref,
            "top1_bitscore": top1_bitscore,
            "top2_bitscore": top2_bitscore,
            "gap_ratio": gap_ratio,
            "top1_matches_expected": top1_matches_expected,
            "support_level": support_level,
        })

        summary_lines.append(f"Candidate: {candidate_id}")
        summary_lines.append(f"  Expected query: {expected_query_id}")
        summary_lines.append(f"  Reciprocal coverage cutoffs: qcov>={rqcov_cutoff}, scov>={rscov_cutoff}")
        summary_lines.append(f"  Total reciprocal hits: {len(sub)}")
        summary_lines.append(f"  Valid reciprocal hits after coverage filter: {len(filtered)}")
        summary_lines.append(f"  Top1 ref: {top1_ref} | bitscore={top1_bitscore}")
        summary_lines.append(f"  Top2 ref: {top2_ref} | bitscore={top2_bitscore}")
        summary_lines.append(f"  Gap ratio: {gap_ratio}")
        summary_lines.append(f"  Top1 matches expected: {top1_matches_expected}")
        summary_lines.append(f"  Support level: {support_level}")
        summary_lines.append("")

    support_df = pd.DataFrame(rows)
    support_df.to_csv(support_file, sep="\t", index=False)

    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    strong_count = int((support_df["support_level"] == "strong").sum())
    ambiguous_count = int((support_df["support_level"] == "ambiguous").sum())
    no_valid_hit_count = int((support_df["support_level"] == "no_valid_hit").sum())

    return {
        "candidates_evaluated": len(support_df),
        "strong_count": strong_count,
        "ambiguous_count": ambiguous_count,
        "no_valid_hit_count": no_valid_hit_count,
        "support_file": str(support_file),
        "summary_file": str(summary_file),
    }

