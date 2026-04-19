from pathlib import Path
import pandas as pd


def identify_first_tier_candidates(config: dict) -> dict:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    top_hits_file = output_dir / "forward_top_hits.tsv"
    candidates_file = output_dir / "forward_candidates.tsv"
    summary_file = output_dir / "forward_breakpoint_summary.txt"

    qcov_cutoff = float(config["thresholds"]["query_cov_cutoff"])
    scov_cutoff = float(config["thresholds"]["subject_cov_cutoff"])
    drop_cutoff = float(config["thresholds"]["breakpoint_drop_cutoff"])

    if not top_hits_file.exists():
        raise FileNotFoundError(f"Top hits file not found: {top_hits_file}")

    df = pd.read_csv(top_hits_file, sep="\t")

    if df.empty:
        empty_cols = list(df.columns) + ["passes_cov_filter", "score_drop_to_next", "selected_as_candidate"]
        pd.DataFrame(columns=empty_cols).to_csv(candidates_file, sep="\t", index=False)
        with open(summary_file, "w", encoding="utf-8") as f:
            f.write("No hits available after forward BLAST preprocessing.\n")
        return {
            "queries_processed": 0,
            "queries_with_candidates": 0,
            "candidates_file": str(candidates_file),
            "summary_file": str(summary_file),
        }

    df["passes_cov_filter"] = (df["qcov"] >= qcov_cutoff) & (df["scov"] >= scov_cutoff)

    selected_blocks = []
    summary_lines = []

    for qseqid, sub in df.groupby("qseqid", sort=False):
        sub = sub.sort_values(
            by=["bitscore", "evalue", "qcov", "scov"],
            ascending=[False, True, False, False]
        ).copy()

        filtered = sub[sub["passes_cov_filter"]].copy()

        summary_lines.append(f"Query: {qseqid}")
        summary_lines.append(f"  Total hits before coverage filter: {len(sub)}")
        summary_lines.append(f"  Hits after coverage filter: {len(filtered)}")
        summary_lines.append(f"  Coverage cutoffs: qcov>={qcov_cutoff}, scov>={scov_cutoff}")

        if filtered.empty:
            sub["score_drop_to_next"] = pd.NA
            sub["selected_as_candidate"] = False
            selected_blocks.append(sub)
            summary_lines.append("  No candidates retained: no hits passed coverage filtering.")
            summary_lines.append("")
            continue

        filtered = filtered.reset_index(drop=True)
        drops = []

        for i in range(len(filtered)):
            if i < len(filtered) - 1:
                current_score = filtered.loc[i, "bitscore"]
                next_score = filtered.loc[i + 1, "bitscore"]
                if current_score == 0:
                    drop = pd.NA
                else:
                    drop = (current_score - next_score) / current_score
            else:
                drop = pd.NA
            drops.append(drop)

        filtered["score_drop_to_next"] = drops

        breakpoint_index = None
        for i, drop in enumerate(drops[:-1]):
            if pd.notna(drop) and drop >= drop_cutoff:
                breakpoint_index = i
                break

        if breakpoint_index is not None:
            candidate_count = breakpoint_index + 1
            summary_lines.append(f"  Breakpoint detected after rank {candidate_count}")
            summary_lines.append(f"  Drop threshold: {drop_cutoff}")
            summary_lines.append(f"  Candidate count retained: {candidate_count}")
        else:
            candidate_count = len(filtered)
            summary_lines.append("  No clear breakpoint detected.")
            summary_lines.append(f"  All coverage-passing hits retained: {candidate_count}")

        filtered["selected_as_candidate"] = False
        filtered.loc[:candidate_count - 1, "selected_as_candidate"] = True

        filtered["expected_query_id"] = qseqid

        sub = sub.merge(
            filtered[["sseqid", "expected_query_id", "score_drop_to_next", "selected_as_candidate"]],
            on="sseqid",
            how="left"
        )

        sub["selected_as_candidate"] = sub["selected_as_candidate"].fillna(False)
        sub["expected_query_id"] = sub["expected_query_id"].fillna(qseqid)
        selected_blocks.append(sub)

        summary_lines.append("  Retained candidates:")
        for _, row in filtered[filtered["selected_as_candidate"]].iterrows():
            summary_lines.append(
                f"    - {row['sseqid']} | bitscore={row['bitscore']} | qcov={row['qcov']:.3f} | scov={row['scov']:.3f}"
            )
        summary_lines.append("")

    final_df = pd.concat(selected_blocks, ignore_index=True)
    final_df.to_csv(candidates_file, sep="\t", index=False)

    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    queries_processed = df["qseqid"].nunique()
    queries_with_candidates = final_df.groupby("qseqid")["selected_as_candidate"].any().sum()

    return {
        "queries_processed": int(queries_processed),
        "queries_with_candidates": int(queries_with_candidates),
        "candidates_file": str(candidates_file),
        "summary_file": str(summary_file),
    }
