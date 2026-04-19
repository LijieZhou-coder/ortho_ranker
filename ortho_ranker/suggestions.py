from pathlib import Path


def build_makeblastdb_command(input_fasta: str, db_prefix: str, title: str) -> str:
    fasta = Path(input_fasta).expanduser()
    prefix = Path(db_prefix).expanduser()

    return (
        f'makeblastdb '
        f'-in "{fasta}" '
        f'-dbtype prot '
        f'-out "{prefix}" '
        f'-title "{title}"'
    )


def generate_blastdb_suggestions(config: dict, blastdb_status: dict) -> list:
    suggestions = []

    target_proteome = config["paths"]["target_proteome"]
    target_prefix = config["paths"]["target_blastdb_prefix"]
    reference_proteome = config["paths"]["reference_proteome"]
    reference_prefix = config["paths"]["reference_blastdb_prefix"]

    project_name = config["project"]["name"]

    if not blastdb_status["target_blastdb"]["all_present"]:
        suggestions.append({
            "name": "build_target_blastdb",
            "command": build_makeblastdb_command(
                target_proteome,
                target_prefix,
                f"{project_name}_target_protein_db"
            )
        })

    if not blastdb_status["reference_blastdb"]["all_present"]:
        suggestions.append({
            "name": "build_reference_blastdb",
            "command": build_makeblastdb_command(
                reference_proteome,
                reference_prefix,
                f"{project_name}_reference_protein_db"
            )
        })

    return suggestions
