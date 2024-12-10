#!/bin/bash

# Check if a directory is provided, otherwise use the current directory
DIRECTORY="${1:-.}"

# Use find to locate all PNG files and convert them to JPG
find "$DIRECTORY" -type f -name '*.png' | while read -r file; do
    # Create the output filename by changing the extension
    output="${file%.png}.jpg"

    # Convert PNG to JPG using ImageMagick's convert command
    magick "$file" "$output"

    # Optional: Print a message indicating success
    echo "Converted '$file' to '$output'"
done
