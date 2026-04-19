from pathlib import Path
from datetime import datetime

from ortho_ranker.config import load_config, validate_config, get_path_status
from ortho_ranker.blastdb import check_blastdb
from ortho_ranker.suggestions import generate_blastdb_suggestions


def run_assessment(config_path: Path) -> dict:
    config = load_config(config_path)
    config_errors = validate_config(config)

    output_dir = Path(config["paths"]["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    path_status = get_path_status(config)

    blastdb_status = {
        "target_blastdb": check_blastdb(config["paths"]["target_blastdb_prefix"]),
        "reference_blastdb": check_blastdb(config["paths"]["reference_blastdb_prefix"]),
    }

    suggestions = generate_blastdb_suggestions(config, blastdb_status)

    assessment = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "config_path": str(config_path),
        "project_name": config["project"]["name"] if "project" in config and "name" in config["project"] else "UNKNOWN",
        "config_errors": config_errors,
        "path_status": path_status,
        "blastdb_status": blastdb_status,
        "suggestions": suggestions,
    }

    write_text_report(output_dir, assessment)
    write_tsv_report(output_dir, assessment)
    write_shell_script(output_dir, assessment)

    return assessment


def write_text_report(output_dir: Path, assessment: dict) -> None:
    report_path = output_dir / "assessment_summary.txt"

    lines = []
    lines.append(f"Assessment time: {assessment['timestamp']}")
    lines.append(f"Config path: {assessment['config_path']}")
    lines.append(f"Project name: {assessment['project_name']}")
    lines.append("")

    if assessment["config_errors"]:
        lines.append("Config validation: FAILED")
        for err in assessment["config_errors"]:
            lines.append(f"- {err}")
    else:
        lines.append("Config validation: PASSED")

    lines.append("")
    lines.append("Path status:")

    for key, info in assessment["path_status"].items():
        if key == "output_dir":
            lines.append(
                f"- {key}: {info['path']} | exists={info['exists']} | is_dir={info['is_dir']}"
            )
        else:
            lines.append(
                f"- {key}: {info['path']} | exists={info['exists']} | is_file={info['is_file']}"
            )

    lines.append("")
    lines.append("BLAST DB status:")

    for db_name, db_info in assessment["blastdb_status"].items():
        lines.append(f"- {db_name}: prefix={db_info['prefix']} | all_present={db_info['all_present']}")
        for suffix, file_info in db_info["files"].items():
            lines.append(
                f"  - {suffix}: {file_info['path']} | exists={file_info['exists']} | is_file={file_info['is_file']}"
            )

    lines.append("")
    lines.append("Suggested commands:")

    if assessment["suggestions"]:
        for item in assessment["suggestions"]:
            lines.append(f"- {item['name']}:")
            lines.append(f"  {item['command']}")
    else:
        lines.append("- No BLAST DB build commands needed.")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def write_tsv_report(output_dir: Path, assessment: dict) -> None:
    report_path = output_dir / "assessment_summary.tsv"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("section\titem\tpath\texists\ttype_ok\textra\n")

        for key, info in assessment["path_status"].items():
            if key == "output_dir":
                type_ok = info["is_dir"] if info["exists"] else True
                f.write(f"path\t{key}\t{info['path']}\t{info['exists']}\t{type_ok}\t.\n")
            else:
                f.write(f"path\t{key}\t{info['path']}\t{info['exists']}\t{info['is_file']}\t.\n")

        for db_name, db_info in assessment["blastdb_status"].items():
            f.write(f"blastdb\t{db_name}\t{db_info['prefix']}\t{db_info['all_present']}\t.\tprefix\n")
            for suffix, file_info in db_info["files"].items():
                f.write(
                    f"blastdb\t{db_name}.{suffix}\t{file_info['path']}\t{file_info['exists']}\t{file_info['is_file']}\tfile\n"
                )

        for item in assessment["suggestions"]:
            f.write(f"suggestion\t{item['name']}\t.\t.\t.\t{item['command']}\n")


def write_shell_script(output_dir: Path, assessment: dict) -> None:
    script_path = output_dir / "build_blastdb_commands.sh"

    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]

    if assessment["suggestions"]:
        for item in assessment["suggestions"]:
            lines.append(item["command"])
    else:
        lines.append('echo "No BLAST DB build commands needed."')

    with open(script_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    script_path.chmod(0o755)
