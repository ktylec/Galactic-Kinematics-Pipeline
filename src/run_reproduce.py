"""Convenience entrypoint to run the original UHST v15 reproduction script.

This preserves the upstream `reproduce_v15.py` as-is and calls it from the repo root.
"""

import runpy
from pathlib import Path

if __name__ == "__main__":
    # Ensure working directory is repo root
    root = Path(__file__).resolve().parents[1]
    # Execute the upstream script
    runpy.run_path(str(root / "src" / "reproduce_v15.py"), run_name="__main__")
