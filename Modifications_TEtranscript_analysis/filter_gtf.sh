#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 input_gtf_file"
    exit 1
fi

# Input file
input_file="$1"

# Extract rows with "transcript" in the 3rd column
awk '$3 == "transcript"' "$input_file"
