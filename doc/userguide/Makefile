# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_userguide.md 010_gettingstarted.md 020_codeoptions.md 021_parameter.md  030_workflow.md
ifneq ($(BUILDALL),)
  SRC += 040_codedocumentation.md
endif
SRC +=  050_tutorials.md tutorials/*/*.md 060_regression_check.md 070_unittest.md 080_installationguidelines.md 081_advanced_installation.md
ifneq ($(BUILDALL),)
  SRC += missing.md
endif
  SRC += references.md

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
	pandoc $(SRC)   -o userguide.pdf -s -N --filter pandoc-citeproc --from markdown+implicit_figures --bibliography userguide.bib --listings --template=template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N --filter pandoc-citeproc --from markdown+implicit_figures --bibliography userguide.bib --listings --template=template.tex

.PHONY: userguide

