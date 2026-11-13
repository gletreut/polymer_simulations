#!/usr/bin/env bash
# Golden-run regression: run bin/simu_langevin on config.yaml in a temporary
# directory and compare thermo.dat and positions_final.dat with the stored
# copies (tolerance 5e-8 on |diff|/max(|ref|, 1)). Exit nonzero on mismatch.
#   compare.sh            compare against stored files
#   compare.sh --update   regenerate the stored files
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
bin="$here/../../bin/simu_langevin"
files="thermo.dat positions_final.dat"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
(cd "$tmp" && "$bin" < "$here/config.yaml" > simu.log 2>&1) || {
  tail -n 20 "$tmp/simu.log"; echo "golden: simulation failed"; exit 1; }

if [[ "${1:-}" == "--update" ]]; then
  for f in $files; do cp "$tmp/$f" "$here/$f"; done
  echo "golden: stored files updated"
  exit 0
fi

python3 - "$here" "$tmp" $files <<'EOF'
import sys
import numpy as np

ref_dir, new_dir, files = sys.argv[1], sys.argv[2], sys.argv[3:]
ok = True
for f in files:
    ref = np.loadtxt(f"{ref_dir}/{f}", ndmin=2)
    new = np.loadtxt(f"{new_dir}/{f}", ndmin=2)
    if ref.shape != new.shape:
        print(f"{f:22s} FAIL  shape {new.shape} != {ref.shape}")
        ok = False
        continue
    # values are printed with about 8 digits: allow a few units in the last one
    err = (np.abs(new - ref) / np.maximum(np.abs(ref), 1.0)).max()
    status = "PASS" if err <= 5e-8 else "FAIL"
    ok &= status == "PASS"
    print(f"{f:22s} max diff {err:.3e}  {status}")
sys.exit(0 if ok else 1)
EOF
