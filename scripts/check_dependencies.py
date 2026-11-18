import importlib.util
import shutil
import sys

PYTHON_PACKAGES = ["decoupler", "omnipath", "pandas", "numpy", "sklearn"]
OPTIONAL_PACKAGES = ["sksurv", "GEOparse", "snakemake"]
OPTIONAL_COMMANDS = ["java", "gdc-client"]


def available_package(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def main() -> int:
    missing_required = [pkg for pkg in PYTHON_PACKAGES if not available_package(pkg)]
    if missing_required:
        print(f"Missing required Python packages: {', '.join(missing_required)}", file=sys.stderr)
        return 1
    for pkg in OPTIONAL_PACKAGES:
        print(f"package\t{pkg}\t{available_package(pkg)}")
    for cmd in OPTIONAL_COMMANDS:
        print(f"command\t{cmd}\t{shutil.which(cmd) is not None}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
