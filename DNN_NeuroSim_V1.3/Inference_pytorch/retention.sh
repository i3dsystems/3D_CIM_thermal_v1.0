#!/bin/bash

time=(1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 3.154e8 1e9)
#drift_coefficient=(1.5345e-14 1.931818e-5 8.6459e-6)
drift_coefficient=(1.5345e-14)

for v in ${drift_coefficient[@]}; do
for t in ${time[@]}; do
    python inference.py --inference 1 --cellBit 1 --onoffratio 100 --t 1e1 --v 1.5345e-14 --detect 1 --target -1
done
done