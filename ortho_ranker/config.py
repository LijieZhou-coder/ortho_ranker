from pathlib import Path
import yaml


def load_config(config_path: Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict):
        raise ValueError("Config file is empty or not a valid YAML mapping.")

    return data


def validate_required_keys(config: dict) -> list:
    errors = []

    required_top_keys = ["project", "paths", "thresholds", "domain", "family", "expression"]
    for key in required_top_keys:
        if key not in config:
            errors.append(f"Missing top-level key: {key}")

    if errors:
        return errors

    required_project_keys = [
        "name",
        "mode",
        "goal",
        "allow_no_unique_hit",
        "return_single_best_only_if_unique",
    ]
    for key in required_project_keys:
        if key not in config["project"]:
            errors.append(f"Missing project.{key}")

    required_paths_keys = [
        "query_fasta",
        "target_proteome",
        "target_blastdb_prefix",
        "reference_proteome",
        "reference_blastdb_prefix",
        "pfam_db",
        "output_dir",
    ]
    for key in required_paths_keys:
        if key not in config["paths"]:
            errors.append(f"Missing paths.{key}")

    required_threshold_keys = [
        "evalue_cutoff",
        "query_cov_cutoff",
        "subject_cov_cutoff",
        "reciprocal_query_cov_cutoff",
        "reciprocal_subject_cov_cutoff",
        "top1_ratio_cutoff",
        "breakpoint_drop_cutoff",
        "breakpoint_median_multiplier",
        "reciprocal_top_gap_cutoff",
        "max_hits_to_keep",
    ]
    for key in required_threshold_keys:
        if key not in config["thresholds"]:
            errors.append(f"Missing thresholds.{key}")

    required_domain_keys = [
        "require_core_domain",
        "require_domain_architecture_match",
        "allow_domain_partial_for_single_best",
    ]
    for key in required_domain_keys:
        if key not in config["domain"]:
            errors.append(f"Missing domain.{key}")

    required_family_keys = ["mode", "rationale"]
    for key in required_family_keys:
        if key not in config["family"]:
            errors.append(f"Missing family.{key}")
    required_expression_keys = [
        "enabled",
        "matrix_file",
        "gene_id_column",
        "expression_metric",
        "target_samples",
        "background_samples",
        "ranking_strategy",
        "min_expression_threshold",
        "include_no_valid_hit",
        "pseudocount",
    ]
    for key in required_expression_keys:
        if key not in config["expression"]:
            errors.append(f"Missing expression.{key}")

    return errors


def _validate_string(config: dict, section: str, key: str, errors: list) -> None:
    value = config[section][key]
    if not isinstance(value, str) or not value.strip():
        errors.append(
            f"{section}.{key} must be a non-empty string, got {type(value).__name__}: {value}"
        )


def _validate_nonempty_list_of_strings(config: dict, section: str, key: str, errors: list) -> None:
    value = config[section][key]
    if not isinstance(value, list) or len(value) == 0:
        errors.append(f"{section}.{key} must be a non-empty list of strings, got {type(value).__name__}: {value}")
        return

    for item in value:
        if not isinstance(item, str) or not item.strip():
            errors.append(f"{section}.{key} must contain only non-empty strings, got item: {item}")
            break


def _validate_bool(config: dict, section: str, key: str, errors: list) -> None:
    value = config[section][key]
    if not isinstance(value, bool):
        errors.append(
            f"{section}.{key} must be a boolean (true/false), got {type(value).__name__}: {value}"
        )


def _validate_number(config: dict, section: str, key: str, errors: list) -> None:
    value = config[section][key]
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        errors.append(
            f"{section}.{key} must be a number, got {type(value).__name__}: {value}"
        )


def _validate_int(config: dict, section: str, key: str, errors: list) -> None:
    value = config[section][key]
    if not isinstance(value, int) or isinstance(value, bool):
        errors.append(
            f"{section}.{key} must be an integer, got {type(value).__name__}: {value}"
        )


def _validate_range(
    config: dict,
    section: str,
    key: str,
    min_value: float | None,
    max_value: float | None,
    errors: list,
) -> None:
    value = config[section][key]
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return

    if min_value is not None and value < min_value:
        errors.append(f"{section}.{key} must be >= {min_value}, got {value}")

    if max_value is not None and value > max_value:
        errors.append(f"{section}.{key} must be <= {max_value}, got {value}")


def _validate_allowed_values(
    config: dict,
    section: str,
    key: str,
    allowed_values: set,
    errors: list,
) -> None:
    value = config[section][key]
    if value not in allowed_values:
        allowed_text = ", ".join(sorted(str(v) for v in allowed_values))
        errors.append(
            f"{section}.{key} must be one of: {allowed_text}; got {value}"
        )


def validate_semantics(config: dict) -> list:
    errors = []

    bool_fields = [
        ("project", "allow_no_unique_hit"),
        ("project", "return_single_best_only_if_unique"),
        ("domain", "require_core_domain"),
        ("domain", "require_domain_architecture_match"),
        ("domain", "allow_domain_partial_for_single_best"),
    ]
    for section, key in bool_fields:
        _validate_bool(config, section, key, errors)

    numeric_fields = [
        ("thresholds", "evalue_cutoff"),
        ("thresholds", "query_cov_cutoff"),
        ("thresholds", "subject_cov_cutoff"),
        ("thresholds", "reciprocal_query_cov_cutoff"),
        ("thresholds", "reciprocal_subject_cov_cutoff"),
        ("thresholds", "top1_ratio_cutoff"),
        ("thresholds", "breakpoint_drop_cutoff"),
        ("thresholds", "breakpoint_median_multiplier"),
        ("thresholds", "reciprocal_top_gap_cutoff"),
    ]
    for section, key in numeric_fields:
        _validate_number(config, section, key, errors)

    _validate_int(config, "thresholds", "max_hits_to_keep", errors)

    zero_to_one_fields = [
        ("thresholds", "query_cov_cutoff"),
        ("thresholds", "subject_cov_cutoff"),
        ("thresholds", "reciprocal_query_cov_cutoff"),
        ("thresholds", "reciprocal_subject_cov_cutoff"),
        ("thresholds", "top1_ratio_cutoff"),
        ("thresholds", "breakpoint_drop_cutoff"),
        ("thresholds", "reciprocal_top_gap_cutoff"),
    ]
    for section, key in zero_to_one_fields:
        _validate_range(config, section, key, 0.0, 1.0, errors)

    _validate_range(config, "thresholds", "evalue_cutoff", 0.0, None, errors)
    _validate_range(config, "thresholds", "breakpoint_median_multiplier", 0.0, None, errors)
    _validate_range(config, "thresholds", "max_hits_to_keep", 1, None, errors)

    allowed_project_modes = {"strict"}
    allowed_project_goals = {"high_precision"}
    allowed_family_modes = {"simple"}

    _validate_allowed_values(config, "project", "mode", allowed_project_modes, errors)
    _validate_allowed_values(config, "project", "goal", allowed_project_goals, errors)
    _validate_allowed_values(config, "family", "mode", allowed_family_modes, errors)

    _validate_bool(config, "expression", "enabled", errors)
    _validate_string(config, "expression", "matrix_file", errors)
    _validate_string(config, "expression", "gene_id_column", errors)
    _validate_string(config, "expression", "expression_metric", errors)
    _validate_string(config, "expression", "ranking_strategy", errors)
    _validate_nonempty_list_of_strings(config, "expression", "target_samples", errors)
    _validate_nonempty_list_of_strings(config, "expression", "background_samples", errors)
    _validate_number(config, "expression", "min_expression_threshold", errors)
    _validate_bool(config, "expression", "include_no_valid_hit", errors)
    _validate_number(config, "expression", "pseudocount", errors)

    _validate_range(config, "expression", "min_expression_threshold", 0.0, None, errors)
    _validate_range(config, "expression", "pseudocount", 0.0, None, errors)

    allowed_expression_metrics = {"FPKM"}
    allowed_ranking_strategies = {"combined"}

    _validate_allowed_values(
        config, "expression", "expression_metric", allowed_expression_metrics, errors
    )
    _validate_allowed_values(
        config, "expression", "ranking_strategy", allowed_ranking_strategies, errors
    )

    return errors


def validate_config(config: dict) -> list:
    required_key_errors = validate_required_keys(config)
    if required_key_errors:
        return required_key_errors

    semantic_errors = validate_semantics(config)
    return semantic_errors


def get_path_status(config: dict) -> dict:
    paths = config["paths"]

    status = {}
    for key in ["query_fasta", "target_proteome", "reference_proteome", "pfam_db"]:
        p = Path(paths[key]).expanduser()
        status[key] = {
            "path": str(p),
            "exists": p.exists(),
            "is_file": p.is_file(),
        }

    output_dir = Path(paths["output_dir"]).expanduser()
    status["output_dir"] = {
        "path": str(output_dir),
        "exists": output_dir.exists(),
        "is_dir": output_dir.is_dir() if output_dir.exists() else False,
    }

    return status


def load_and_validate_config(config_path: Path) -> dict:
    config = load_config(config_path)
    errors = validate_config(config)

    if errors:
        error_text = "\n".join(errors)
        raise ValueError(error_text)

    return config
