set -euo pipefail
cd "$(dirname "$0")/.."
make -j"$(sysctl -n hw.ncpu)" OMP=0
make hs-run
