#!/bin/bash

# Go to your Bubbles folder
cd ~/Desktop/Bulles/Bubbles

# Create output folder for .dat files
mkdir -p dat_files

count=0

for file in *.vtk; do
    num=$(printf "%04d" $count)
    outFile="dat_files/dat_$num.dat"
    echo "Converting $file -> $outFile ..."
    
    # Run your converter
    ./conv "$file" "$outFile"
    
    count=$((count + 1))
done    