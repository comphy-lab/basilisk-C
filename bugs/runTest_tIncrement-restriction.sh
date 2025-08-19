#!/bin/bash

echo "=== Testing tfake=1.2e-7 (should work) ==="
qcc -O2 -Wall -disable-dimensions -DTFAKE=1.2e-7 tIncrement-restriction.c -o tIncrement-restriction -lm
./tIncrement-restriction | head -20

echo -e "\n=== Testing tfake=1.19e-7 (will fail) ==="
qcc -O2 -Wall -disable-dimensions -DTFAKE=1.19e-7 tIncrement-restriction.c -o tIncrement-restriction -lm
./tIncrement-restriction 2>&1 | head -20