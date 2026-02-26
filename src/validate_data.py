import argparse
from pathlib import Path

EXPECTED_COLS = 6  # R, Vobs, eV, Vgas, Vdisk, Vbul

def is_data_row(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#"):
        return False
    parts = s.split()
    if len(parts) < EXPECTED_COLS:
        return False
    try:
        [float(x) for x in parts[:EXPECTED_COLS]]
        return True
    except ValueError:
        return False

def main():
    ap = argparse.ArgumentParser(description="Quick sanity-check for SPARC rotmod .dat files.")
    ap.add_argument("--data_dir", default="data/sparc_dat", help="Folder with *_rotmod.dat files.")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    files = sorted(data_dir.glob("*_rotmod.dat"))
    if not files:
        raise SystemExit(f"No *_rotmod.dat files found in {data_dir}")

    bad = []
    total_rows = 0
    for fp in files:
        rows = 0
        with fp.open("r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if is_data_row(ln):
                    rows += 1
        total_rows += rows
        if rows == 0:
            bad.append(fp.name)

    print(f"Files: {len(files)}")
    print(f"Total numeric rows (>= {EXPECTED_COLS} cols): {total_rows}")
    if bad:
        print("Files with zero parsed rows:")
        for b in bad:
            print(" -", b)
        raise SystemExit(2)

if __name__ == "__main__":
    main()
