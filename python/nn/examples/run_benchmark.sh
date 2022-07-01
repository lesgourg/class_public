#!/usr/bin/env sh

# This bash script can be used to run benchmarks on the calculation time.
# It is to be executed in the terminal via the command:
# > bash run_benchmark.sh


# first run benchmark loop
WARMUP=5
ITERATIONS=50
THREAD_COUNTS="1 2 4 8 12 16"

for i in $THREAD_COUNTS
do
    OMP_NUM_THREADS=$i python3 -m benchmark $i -w $WARMUP -i $ITERATIONS
done

# After benchmark run the results can be plotted using a plotting script.
# They are stored in the workspace directory: path/to/workspace/results/benchmark/
python3 ../plotting/plot_benchmark.py