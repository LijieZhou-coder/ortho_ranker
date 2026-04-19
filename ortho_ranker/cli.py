from pathlib import Path
import pandas as pd
import typer

from ortho_ranker.assess import run_assessment
from ortho_ranker.config import load_and_validate_config
from ortho_ranker.runner import run_forward_blast, preprocess_forward_hits
from ortho_ranker.candidates import identify_first_tier_candidates
from ortho_ranker.reciprocal import (
    export_candidate_fasta,
    run_reciprocal_blast,
    summarize_reciprocal_support,
)
from ortho_ranker.expression import integrate_expression_evidence

app = typer.Typer(help="Ortho Ranker: high-precision homolog screening tool")


@app.command()
def assess(
    config: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
):
    result = run_assessment(config)

    typer.echo(f"[ASSESS] project = {result['project_name']}")
    typer.echo(f"[ASSESS] config path = {result['config_path']}")

    if result["config_errors"]:
        typer.echo("[ASSESS] config validation = FAILED")
        for err in result["config_errors"]:
            typer.echo(f"  - {err}")
    else:
        typer.echo("[ASSESS] config validation = PASSED")

    typer.echo("[ASSESS] path status:")
    for key, info in result["path_status"].items():
        if key == "output_dir":
            typer.echo(
                f"  - {key}: exists={info['exists']}, is_dir={info['is_dir']}, path={info['path']}"
            )
        else:
            typer.echo(
                f"  - {key}: exists={info['exists']}, is_file={info['is_file']}, path={info['path']}"
            )

    typer.echo("[ASSESS] BLAST DB status:")
    for db_name, db_info in result["blastdb_status"].items():
        typer.echo(
            f"  - {db_name}: all_present={db_info['all_present']}, prefix={db_info['prefix']}"
        )
        for suffix, file_info in db_info["files"].items():
            typer.echo(
                f"      * {suffix}: exists={file_info['exists']}, is_file={file_info['is_file']}, path={file_info['path']}"
            )

    typer.echo("[ASSESS] suggested commands:")
    if result["suggestions"]:
        for item in result["suggestions"]:
            typer.echo(f"  - {item['name']}: {item['command']}")
    else:
        typer.echo("  - No BLAST DB build commands needed.")

    typer.echo("[ASSESS] reports written to output_dir")
    typer.echo("[ASSESS] shell script written: build_blastdb_commands.sh")


@app.command()
def run(
    config: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
):
    try:
        cfg = load_and_validate_config(config)
    except ValueError as e:
        typer.echo("[RUN] config validation = FAILED")
        for line in str(e).splitlines():
            typer.echo(f"  {line}")
        raise typer.Exit(code=2)

    result = run_forward_blast(cfg)

    typer.echo(f"[RUN] blastp command = {result['command']}")
    typer.echo(f"[RUN] returncode = {result['returncode']}")

    if result["stdout"]:
        typer.echo("[RUN] stdout:")
        typer.echo(result["stdout"])

    if result["stderr"]:
        typer.echo("[RUN] stderr:")
        typer.echo(result["stderr"])

    if result["returncode"] != 0:
        raise typer.Exit(code=result["returncode"])

    typer.echo(f"[RUN] forward BLAST output = {result['output_file']}")

    top_result = preprocess_forward_hits(cfg)
    typer.echo(f"[RUN] raw rows = {top_result['raw_rows']}")
    typer.echo(f"[RUN] top rows = {top_result['top_rows']}")
    typer.echo(f"[RUN] top hits file = {top_result['top_hits_file']}")

    candidate_result = identify_first_tier_candidates(cfg)
    typer.echo(f"[RUN] queries processed = {candidate_result['queries_processed']}")
    typer.echo(f"[RUN] queries with candidates = {candidate_result['queries_with_candidates']}")
    typer.echo(f"[RUN] candidate file = {candidate_result['candidates_file']}")
    typer.echo(f"[RUN] breakpoint summary = {candidate_result['summary_file']}")

    export_result = export_candidate_fasta(cfg)
    typer.echo(f"[RUN] reciprocal candidate count = {export_result['candidate_count']}")
    typer.echo(f"[RUN] reciprocal exported count = {export_result['exported_count']}")
    typer.echo(f"[RUN] reciprocal fasta = {export_result['reciprocal_fasta']}")

    reciprocal_result = run_reciprocal_blast(cfg)
    typer.echo(f"[RUN] reciprocal blastp command = {reciprocal_result['command']}")
    typer.echo(f"[RUN] reciprocal returncode = {reciprocal_result['returncode']}")

    if reciprocal_result["stdout"]:
        typer.echo("[RUN] reciprocal stdout:")
        typer.echo(reciprocal_result["stdout"])

    if reciprocal_result["stderr"]:
        typer.echo("[RUN] reciprocal stderr:")
        typer.echo(reciprocal_result["stderr"])

    if reciprocal_result["returncode"] != 0:
        raise typer.Exit(code=reciprocal_result["returncode"])

    typer.echo(f"[RUN] reciprocal BLAST output = {reciprocal_result['output_file']}")

    support_result = summarize_reciprocal_support(cfg)
    typer.echo(f"[RUN] reciprocal candidates evaluated = {support_result['candidates_evaluated']}")
    typer.echo(f"[RUN] reciprocal strong count = {support_result['strong_count']}")
    typer.echo(f"[RUN] reciprocal ambiguous count = {support_result['ambiguous_count']}")
    typer.echo(f"[RUN] reciprocal no_valid_hit count = {support_result.get('no_valid_hit_count', 0)}")
    typer.echo(f"[RUN] reciprocal support file = {support_result['support_file']}")
    typer.echo(f"[RUN] reciprocal summary file = {support_result['summary_file']}")

    if cfg["expression"]["enabled"]:
        support_df = pd.read_csv(support_result["support_file"], sep="\t")
        expr_result = integrate_expression_evidence(support_df, cfg)

        typer.echo("[RUN] expression ranking enabled = True")
        typer.echo(f"[RUN] expression ranked candidates = {expr_result['ranked_candidates']}")
        typer.echo(f"[RUN] expression top candidate = {expr_result['top_candidate']}")
        typer.echo(f"[RUN] expression ranking file = {expr_result['ranking_file']}")
        typer.echo(f"[RUN] expression summary file = {expr_result['summary_file']}")
    else:
        typer.echo("[RUN] expression ranking enabled = False")


def main():
    app()


if __name__ == "__main__":
    main()
