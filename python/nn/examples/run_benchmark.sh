#!/usr/bin/env sh

WARMUP=2
ITERATIONS=25
THREAD_COUNTS="1 2 4 8 16"

for i in $THREAD_COUNTS
do
    OMP_NUM_THREADS=$i python3 -m nn.examples.benchmark $i -w $WARMUP -i $ITERATIONS
done
