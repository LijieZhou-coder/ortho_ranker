from pathlib import Path
import subprocess

import pandas as pd
import pytest

from ortho_ranker.expression import integrate_expression_evidence


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run_cli(args):
    cmd = ["ortho-ranker"] + args
    return subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )


def test_run_with_expression_outputs_files():
    result = run_cli(["run", "config.test.yaml"])

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    ranking_file = PROJECT_ROOT / "results" / "test_run" / "expression_ranking.tsv"
    summary_file = PROJECT_ROOT / "results" / "test_run" / "expression_summary.txt"

    assert ranking_file.exists(), "expression_ranking.tsv was not created"
    assert summary_file.exists(), "expression_summary.txt was not created"

    ranking_df = pd.read_csv(ranking_file, sep="\t")
    assert "expression_priority_score" in ranking_df.columns
    assert "recommended_by_expression" in ranking_df.columns
    assert "expression_data_found" in ranking_df.columns


def test_expression_top_candidate_and_explainability_summary():
    result = run_cli(["run", "config.test.yaml"])

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    summary_file = PROJECT_ROOT / "results" / "test_run" / "expression_summary.txt"
    content = summary_file.read_text(encoding="utf-8")

    assert "expression_metric: FPKM" in content
    assert "ranking_strategy: combined" in content

    assert "top_candidate: target_hit_1" in content
    assert "top_support_level: ambiguous" in content
    assert "top_target_mean_fpkm: 105.0" in content
    assert "top_background_mean_fpkm: 4.5" in content
    assert "top_specificity_ratio:" in content
    assert "top_expression_priority_score:" in content


def test_expression_ranking_matches_expected_biological_order():
    result = run_cli(["run", "config.test.yaml"])

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    ranking_file = PROJECT_ROOT / "results" / "test_run" / "expression_ranking.tsv"
    df = pd.read_csv(ranking_file, sep="\t")

    assert len(df) >= 2

    top_row = df.iloc[0]
    second_row = df.iloc[1]

    assert top_row["candidate_id"] == "target_hit_1"
    assert int(top_row["expression_rank"]) == 1
    assert bool(top_row["recommended_by_expression"]) is True

    assert second_row["candidate_id"] == "target_hit_2"
    assert int(second_row["expression_rank"]) == 2
    assert bool(second_row["recommended_by_expression"]) is False

    assert top_row["target_mean_fpkm"] > second_row["target_mean_fpkm"]
    assert top_row["specificity_ratio"] > second_row["specificity_ratio"]
    assert top_row["expression_priority_score"] > second_row["expression_priority_score"]


def test_missing_expression_id_candidate_is_retained_and_marked():
    result = run_cli(["run", "config.missing_expression_id.yaml"])

    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    ranking_file = PROJECT_ROOT / "results" / "missing_expression_id_run" / "expression_ranking.tsv"
    summary_file = PROJECT_ROOT / "results" / "missing_expression_id_run" / "expression_summary.txt"

    assert ranking_file.exists()
    assert summary_file.exists()

    df = pd.read_csv(ranking_file, sep="\t")

    missing_row = df.loc[df["candidate_id"] == "target_hit_2"].iloc[0]

    assert bool(missing_row["expression_data_found"]) is False
    assert pd.isna(missing_row["expression_priority_score"])

    summary_text = summary_file.read_text(encoding="utf-8")
    assert "expression_data_missing_count: 1" in summary_text
    assert "missing_expression_candidates: target_hit_2" in summary_text


def test_bad_expression_matrix_missing_required_column_fails_cleanly():
    result = run_cli(["run", "config.bad_expression_columns.yaml"])

    assert result.returncode != 0

    combined_output = result.stdout + "\n" + result.stderr
    assert "Missing expression matrix columns" in combined_output
    assert "petal_2" in combined_output


def test_combined_score_penalizes_low_expression_pseudo_specific_candidate(tmp_path):
    output_dir = tmp_path / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    expr_matrix = tmp_path / "expression.tsv"
    expr_matrix.write_text(
        "\n".join([
            "gene_id\tpetal_1\tpetal_2\tleaf_1\tleaf_2\tGO_annotation\tKEGG",
            "target_hit_1\t100\t110\t5\t4\tchlorophyll biosynthesis\tko00195",
            "target_hit_2\t60\t55\t20\t18\tstress response\tko04075",
            "target_hit_3\t0.2\t0.3\t0.0\t0.0\tpseudo-specific low-expression gene\tko99999",
        ]) + "\n",
        encoding="utf-8",
    )

    cfg = {
        "paths": {
            "output_dir": str(output_dir),
        },
        "expression": {
            "enabled": True,
            "matrix_file": str(expr_matrix),
            "gene_id_column": "gene_id",
            "expression_metric": "FPKM",
            "target_samples": ["petal_1", "petal_2"],
            "background_samples": ["leaf_1", "leaf_2"],
            "ranking_strategy": "combined",
            "min_expression_threshold": 0.0,
            "include_no_valid_hit": True,
            "pseudocount": 1.0,
        },
    }

    candidates_df = pd.DataFrame([
        {"candidate_id": "target_hit_1", "expected_query_id": "query_protein_1", "support_level": "ambiguous"},
        {"candidate_id": "target_hit_2", "expected_query_id": "query_protein_1", "support_level": "ambiguous"},
        {"candidate_id": "target_hit_3", "expected_query_id": "query_protein_1", "support_level": "ambiguous"},
    ])

    result = integrate_expression_evidence(candidates_df, cfg)
    ranked_df = result["ranked_df"]

    assert len(ranked_df) == 3

    target1 = ranked_df.loc[ranked_df["candidate_id"] == "target_hit_1"].iloc[0]
    target2 = ranked_df.loc[ranked_df["candidate_id"] == "target_hit_2"].iloc[0]
    target3 = ranked_df.loc[ranked_df["candidate_id"] == "target_hit_3"].iloc[0]

    assert bool(target1["expression_data_found"]) is True
    assert bool(target2["expression_data_found"]) is True
    assert bool(target3["expression_data_found"]) is True

    assert target1["target_mean_fpkm"] == 105.0
    assert target1["background_mean_fpkm"] == 4.5

    assert target3["target_mean_fpkm"] == 0.25
    assert target3["background_mean_fpkm"] == 0.0

    assert target1["expression_priority_score"] > target2["expression_priority_score"]
    assert target2["expression_priority_score"] > target3["expression_priority_score"]
    assert target1["expression_priority_score"] > target3["expression_priority_score"]

    assert int(target1["expression_rank"]) < int(target3["expression_rank"])


def test_expression_empty_candidates_returns_empty_ranking_and_summary(tmp_path):
    output_dir = tmp_path / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    expr_matrix = tmp_path / "expression.tsv"
    expr_matrix.write_text(
        "\n".join([
            "gene_id\tpetal_1\tpetal_2\tleaf_1\tleaf_2",
            "target_hit_1\t100\t110\t5\t4",
            "target_hit_2\t60\t55\t20\t18",
        ]) + "\n",
        encoding="utf-8",
    )

    cfg = {
        "paths": {
            "output_dir": str(output_dir),
        },
        "expression": {
            "enabled": True,
            "matrix_file": str(expr_matrix),
            "gene_id_column": "gene_id",
            "expression_metric": "FPKM",
            "target_samples": ["petal_1", "petal_2"],
            "background_samples": ["leaf_1", "leaf_2"],
            "ranking_strategy": "combined",
            "min_expression_threshold": 0.0,
            "include_no_valid_hit": True,
            "pseudocount": 1.0,
        },
    }

    empty_candidates_df = pd.DataFrame(
        columns=["candidate_id", "expected_query_id", "support_level"]
    )

    result = integrate_expression_evidence(empty_candidates_df, cfg)

    assert result["ranked_candidates"] == 0
    assert result["top_candidate"] is None

    ranking_file = Path(result["ranking_file"])
    summary_file = Path(result["summary_file"])

    assert ranking_file.exists()
    assert summary_file.exists()

    ranking_df = pd.read_csv(ranking_file, sep="\t")
    assert ranking_df.empty

    summary_text = summary_file.read_text(encoding="utf-8")
    assert "ranked_candidates: 0" in summary_text
    assert "expression_data_found_count: 0" in summary_text
    assert "expression_data_missing_count: 0" in summary_text
    assert "top_candidate: NA" in summary_text
