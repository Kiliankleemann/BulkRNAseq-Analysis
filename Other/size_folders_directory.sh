#!/bin/bash

# Define the output file name
output_file="folder_sizes.txt"

# Clear the output file if it exists
> "$output_file"

# Iterate over each folder in the current directory
for folder in */; do
    # Get the folder name (remove trailing '/')
    folder_name="${folder%/}"

    # Get the size of the folder
    folder_size=$(du -sh "$folder_name" | cut -f1)

    # Write the folder name and size to the output file
    echo "$folder_name: $folder_size" >> "$output_file"
done

echo "Folder sizes have been written to $output_file"
