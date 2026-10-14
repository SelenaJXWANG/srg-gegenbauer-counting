"""Verify the data manifest, then all six exact raw SRG identities."""
from pathlib import Path
from hashlib import sha256
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CASES = ["540_147_66_30", "351_140_73_44", "476_133_60_28",
         "703_182_81_35", "550_162_75_36", "405_132_63_33"]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-files", action="store_true",
                        help="Only check file hashes; do not reprove the identities.")
    args = parser.parse_args()
    checked = 0
    for line in (ROOT / "MANIFEST.sha256").read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = ROOT / name
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit("File hash mismatch: " + name)
        checked += 1
    print("Manifest verified: %d files." % checked, flush=True)
    if args.check_files:
        return
    subprocess.run([sys.executable, str(ROOT / "verify_raw_certificate.py")]
                   + [str(ROOT / (case + ".raw.json")) for case in CASES[:5]], check=True)
    subprocess.run([sys.executable, str(ROOT / "verify_raw_v2_certificate.py"),
                    str(ROOT / (CASES[5] + ".raw.json"))], check=True)
    print("All six exact identities verified.", flush=True)

if __name__ == "__main__":
    main()
