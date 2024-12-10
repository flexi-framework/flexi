#!/bin/bash

# Check prerequisites
python3 -c "import sphinx"

# $? Stores the exit value of the last command that was executed.
if [[ $? -eq 0 ]]; then
  # Compile html files
  python3 -m sphinx -T -E -b html -d _build/doctrees -D language=en . _build/html
  # python3 -m sphinx -T -E -b spelling -d _build/doctrees -D language=en . _build/html

  # Output info where the html file is
  echo -e "\nThe HTML files have been created. Run, e.g., 'firefox _build/html/index.html &' to view the documentation."
else
  echo -e "\nError: Could not build the documentation due to import errors in python! Fix them and run the script again."
fi
