#!/bin/bash

# Check prerequisites
python3 -c "import sphinx"

# $? Stores the exit value of the last command that was executed.
if [[ $? -eq 0 ]]; then
  # Compile latex files
  python3 -m sphinx -b latex -D language=en -d _build/doctrees . _build/latex

  # Switch to latex source files
  cd _build/latex

  # Compile pdf file(s)
  latexmk -r latexmkrc -pdf -f -dvi- -ps- -jobname=flexi -interaction=nonstopmode

  # Output info where the pdf file is
  echo -e "\n The PDF has been created under ./_build/latex/flexi.pdf"
else
  echo -e "\nError: Could not build the documentation due to import errors in python! Fix them and run the script again."
fi
