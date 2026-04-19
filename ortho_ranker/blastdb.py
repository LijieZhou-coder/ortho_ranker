from pathlib import Path


def check_blastdb(prefix: str) -> dict:
    p = Path(prefix).expanduser()

    required_files = {
        "phr": Path(str(p) + ".phr"),
        "pin": Path(str(p) + ".pin"),
        "psq": Path(str(p) + ".psq"),
    }

    file_status = {
        suffix: {
            "path": str(path),
            "exists": path.exists(),
            "is_file": path.is_file(),
        }
        for suffix, path in required_files.items()
    }

    all_present = all(info["exists"] and info["is_file"] for info in file_status.values())

    return {
        "prefix": str(p),
        "all_present": all_present,
        "files": file_status,
    }
