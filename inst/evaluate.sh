#!/bin/bash

set -e -u -o pipefail

cat tryptase.out | awk '{printf "./Tryplik %d %d %d %d %d | sort -nk2 | tail -n 1\n", $6, $7, $8, $9, $10}' > in1

chmod +x in1

./in1 > tryptaseMLE

exit 0
