SHELL := /bin/bash

ENV_NAME := polio-capsid-env
THREADS ?= 4
PYTHON ?= python3
SNAKEMAKE ?= snakemake
CONFIG ?= config/config.yaml

.PHONY: env run dry test validate clean help

help:
	@echo "make env   - create conda env"
	@echo "make run   - run full Snakemake workflow"
	@echo "make dry   - dry run to show planned steps"
	@echo "make test  - quick sanity check on plotting"
	@echo "make validate - compile scripts and run lightweight checks"
	@echo "make clean - remove work and results"

env:
	conda env create -f env/environment.yml || echo "Env may already exist"
	@echo "Activate with: conda activate $(ENV_NAME)"

run:
	@if [ -z "$$NCBI_EMAIL" ]; then echo "Set NCBI_EMAIL before running"; exit 1; fi
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) -c $(THREADS) --printshellcmds

dry:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) -n -c $(THREADS)

test:
	$(PYTHON) -m pip install -r env/requirements.txt
	$(PYTHON) -m compileall -q analysis/scripts
	$(PYTHON) analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
	@echo "Wrote results-example/example_plot.png"

validate: test
	@if command -v $(SNAKEMAKE) >/dev/null 2>&1; then \
		echo "Snakemake version: $$($(SNAKEMAKE) --version)"; \
		echo "Run make dry CONFIG=$(CONFIG) after configuring local FASTQ inputs."; \
	else \
		echo "Snakemake not found; install env/environment.yml before workflow dry-run"; \
	fi

clean:
	rm -rf work results logs refs .snakemake
