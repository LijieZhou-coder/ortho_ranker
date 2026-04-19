from pathlib import Path
import math
import pandas as pd


def load_expression_subset(config: dict) -> pd.DataFrame:
    expr_cfg = config["expression"]
    matrix_file = Path(expr_cfg["matrix_file"]).expanduser()

    df = pd.read_csv(matrix_file, sep=None, engine="python")

    gene_id_col = expr_cfg["gene_id_column"]
    target_samples = expr_cfg["target_samples"]
    background_samples = expr_cfg["background_samples"]

    required_columns = [gene_id_col] + target_samples + background_samples
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        missing_text = ", ".join(missing_columns)
        raise ValueError(f"Missing expression matrix columns: {missing_text}")

    subset = df[required_columns].copy()

    numeric_columns = target_samples + background_samples
    for col in numeric_columns:
        subset[col] = pd.to_numeric(subset[col], errors="coerce")

    return subset


def _combined_score(target_mean: float, background_mean: float, pseudocount: float) -> float | None:
    if pd.isna(target_mean) or pd.isna(background_mean):
        return None

    return math.log2(target_mean + pseudocount) + math.log2(
        (target_mean + pseudocount) / (background_mean + pseudocount)
    )


def rank_candidates_by_expression(
    candidates_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    expr_cfg = config["expression"]
    gene_id_col = expr_cfg["gene_id_column"]
    target_samples = expr_cfg["target_samples"]
    background_samples = expr_cfg["background_samples"]
    pseudocount = expr_cfg["pseudocount"]
    min_expression_threshold = expr_cfg["min_expression_threshold"]
    include_no_valid_hit = expr_cfg["include_no_valid_hit"]

    merged = candidates_df.merge(
        expression_df,
        how="left",
        left_on="candidate_id",
        right_on=gene_id_col,
    )

    if gene_id_col in merged.columns:
        merged = merged.drop(columns=[gene_id_col])

    if not include_no_valid_hit and "support_level" in merged.columns:
        merged = merged[merged["support_level"] != "no_valid_hit"].copy()

    merged["expression_data_found"] = merged[target_samples].notna().any(axis=1)

    merged["target_mean_fpkm"] = merged[target_samples].mean(axis=1, skipna=True)
    merged["background_mean_fpkm"] = merged[background_samples].mean(axis=1, skipna=True)

    merged["specificity_ratio"] = (
        (merged["target_mean_fpkm"] + pseudocount)
        / (merged["background_mean_fpkm"] + pseudocount)
    )

    merged["expression_priority_score"] = merged.apply(
        lambda row: _combined_score(
            row["target_mean_fpkm"],
            row["background_mean_fpkm"],
            pseudocount,
        ),
        axis=1,
    )

    low_expr_mask = (
        merged["expression_data_found"]
        & merged["target_mean_fpkm"].notna()
        & (merged["target_mean_fpkm"] < min_expression_threshold)
    )
    merged.loc[low_expr_mask, "expression_priority_score"] = None

    merged = merged.sort_values(
        by=["expression_data_found", "expression_priority_score"],
        ascending=[False, False],
        na_position="last",
    ).reset_index(drop=True)

    merged["expression_rank"] = range(1, len(merged) + 1)

    merged["recommended_by_expression"] = False
    valid_recommendable = (
        merged["expression_data_found"]
        & merged["expression_priority_score"].notna()
    )
    if valid_recommendable.any():
        top_idx = merged[valid_recommendable]["expression_priority_score"].idxmax()
        merged.loc[top_idx, "recommended_by_expression"] = True

    return merged


def write_expression_outputs(ranked_df: pd.DataFrame, config: dict) -> tuple[str, str]:
    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    ranking_path = output_dir / "expression_ranking.tsv"
    summary_path = output_dir / "expression_summary.txt"

    ranked_df.to_csv(ranking_path, sep="\t", index=False)

    expr_cfg = config["expression"]
    target_samples = expr_cfg["target_samples"]
    background_samples = expr_cfg["background_samples"]

    n_total = len(ranked_df)
    n_found = int(ranked_df["expression_data_found"].sum())
    n_missing = n_total - n_found

    recommended = ranked_df[ranked_df["recommended_by_expression"] == True]

    lines = [
        f"expression_metric: {expr_cfg['expression_metric']}",
        f"ranking_strategy: {expr_cfg['ranking_strategy']}",
        f"target_samples: {', '.join(target_samples)}",
        f"background_samples: {', '.join(background_samples)}",
        f"ranked_candidates: {n_total}",
        f"expression_data_found_count: {n_found}",
        f"expression_data_missing_count: {n_missing}",
    ]

    if not recommended.empty:
        top = recommended.iloc[0]
        lines.extend([
            f"top_candidate: {top['candidate_id']}",
            f"top_support_level: {top.get('support_level', 'NA')}",
            f"top_target_mean_fpkm: {top['target_mean_fpkm']}",
            f"top_background_mean_fpkm: {top['background_mean_fpkm']}",
            f"top_specificity_ratio: {top['specificity_ratio']}",
            f"top_expression_priority_score: {top['expression_priority_score']}",
        ])

        if top.get("support_level", "") == "no_valid_hit":
            lines.append("warning: top candidate is supported by expression but reciprocal support_level is no_valid_hit")
    else:
        lines.append("top_candidate: NA")
        lines.append("warning: no candidate had usable expression data for ranking")

    if n_missing > 0:
        missing_ids = ranked_df.loc[~ranked_df["expression_data_found"], "candidate_id"].tolist()
        lines.append("missing_expression_candidates: " + ", ".join(map(str, missing_ids)))

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return str(ranking_path), str(summary_path)


def integrate_expression_evidence(
    candidates_df: pd.DataFrame,
    config: dict,
) -> dict:
    expression_df = load_expression_subset(config)
    ranked_df = rank_candidates_by_expression(candidates_df, expression_df, config)
    ranking_path, summary_path = write_expression_outputs(ranked_df, config)

    recommended = ranked_df[ranked_df["recommended_by_expression"] == True]
    top_candidate = recommended.iloc[0]["candidate_id"] if not recommended.empty else None

    return {
        "ranked_df": ranked_df,
        "ranked_candidates": len(ranked_df),
        "top_candidate": top_candidate,
        "ranking_file": ranking_path,
        "summary_file": summary_path,
    }
