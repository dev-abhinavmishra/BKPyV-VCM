"""DEPRECATED shim -- use scripts/generate_bkpyv_review_bundle.py.

The previous implementation contained a hardcoded linear copies/mL scaling
and unverifiable cohort-calibration claims. Reviewer-facing outputs now come
from the canonical review bundle (ODE model + assumption-labelled bridge +
manifest with solver provenance), and from scripts/simulate_drug_switch.py
for the switch comparison.
"""

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo = Path(__file__).parent
    print(__doc__)
    print()
    return subprocess.call(
        [sys.executable, str(repo / "scripts" / "generate_bkpyv_review_bundle.py")],
        cwd=repo,
    )


if __name__ == "__main__":
    raise SystemExit(main())
