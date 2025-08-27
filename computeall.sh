#!/bin/bash

cd ~/Desktop/Bulles/Bubbles

# Create result folder if it doesn't exist
mkdir -p result

count=0

for file in *.vtk; do
    num=$(printf "%04d" $count)
    echo "Processing $file -> result/result_$num.dat ..."
    
    # Run MPI program with input and output arguments
    mpirun -np 2 ./keep "$file" "result/result_$num.dat"
    
    count=$((count + 1))
done