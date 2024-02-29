#Renaming files in directory - replace all characters before hyphen with empty

for file in *; do
    # Check if the file name contains a hyphen
    if [[ $file == *-* ]]; then
        # Substitute all characters before and including the hyphen with an empty string
        new_name="${file#*-}"
        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed $file to $new_name"
    else
        echo "Skipping $file (no hyphen found)"
    fi
done
