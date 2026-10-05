#!/usr/bin/env bash

set -u -o pipefail
shopt -s nullglob

cd "$(dirname "${BASH_SOURCE[0]}")"

mkdir -p ../../../../test/data/simulation/fasta/
mkdir -p ../../../../test/data/simulation/ground_truth/

scripts=(create_dna-*.py create_rna-*.py)

for script in "${scripts[@]}"; do
  python ${script}
done
