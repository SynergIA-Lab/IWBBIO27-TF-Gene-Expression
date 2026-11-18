"""Download the curated regulons once (needs internet): CollecTRI and DoRothEA A-C via OmniPath."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from regbio.regulons import load_collectri, load_dorothea


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collectri", required=True)
    parser.add_argument("--dorothea", required=True)
    args = parser.parse_args()
    for path, loader in ((args.collectri, load_collectri), (args.dorothea, load_dorothea)):
        net = loader()
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        net.to_csv(path, sep="\t", index=False)
        print(f"{path}: {net['source'].nunique()} TFs, {len(net)} interactions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
