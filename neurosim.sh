#!/bin/bash

time=(1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 3.154e8 1e9)
drift_coefficient=("$@")

# Path to NeuroSim
cd "DNN_NeuroSim_V1.3/Inference_pytorch"

for v in ${drift_coefficient[@]}; do
for t in ${time[@]}; do
    python inference.py --inference 1 --cellBit 1 --onoffratio 100 --t $t --v $v --detect 1 --target -1
done
done