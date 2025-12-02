# Benchmark harness for VCFX
ROOT_DIR := $(abspath ..)
DATA_DIR := $(ROOT_DIR)/benchmarks/data
RESULTS_DIR := $(ROOT_DIR)/benchmarks/results
BENCH_SCRIPTS := $(ROOT_DIR)/benchmarks/scripts
MK_TASKS := $(ROOT_DIR)/benchmarks/mk_tasks.py
TASKS_FILE := $(ROOT_DIR)/benchmarks/tasks.yaml
RESULTS_CSV := $(ROOT_DIR)/benchmarks/results.csv

VCFX_BIN_DIR ?= $(ROOT_DIR)/build/src
VCFX_VARIANT_COUNTER := $(VCFX_BIN_DIR)/VCFX_variant_counter/VCFX_variant_counter
VCFX_ALLELE_FREQ_CALC := $(VCFX_BIN_DIR)/VCFX_allele_freq_calc/VCFX_allele_freq_calc
VCFX_VARIANT_CLASSIFIER := $(VCFX_BIN_DIR)/VCFX_variant_classifier/VCFX_variant_classifier
VCFX_MISSING_DETECTOR := $(VCFX_BIN_DIR)/VCFX_missing_detector/VCFX_missing_detector
VCFX_RECORD_FILTER := $(VCFX_BIN_DIR)/VCFX_record_filter/VCFX_record_filter
VCFX_VALIDATOR := $(VCFX_BIN_DIR)/VCFX_validator/VCFX_validator

# Comparison tools (if available)
BCFTOOLS ?= bcftools
VCFTOOLS ?= vcftools
GATK ?= gatk
VCFLIB_vcfstats ?= vcfstats
VCFLIB_vcffilter ?= vcffilter

# Real 1000 Genomes Phase 3 data URLs
# Note: SHA checksums not verified - these are large public datasets
URL_chr20 := https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr20.1kg.phase3.v5a.vcf.gz
URL_chr21 := https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr21.1kg.phase3.v5a.vcf.gz
URL_intervals := https://raw.githubusercontent.com/samtools/hts-specs/master/test/tabix/example.bed

.PHONY: all datasets tools-check install-deps clean

all: datasets tools-check $(RESULTS_CSV)

$(RESULTS_CSV): harness.mk tasks.mk
	$(MAKE) -f harness.mk RESULTS_CSV=$(RESULTS_CSV) BENCH_SCRIPTS=$(BENCH_SCRIPTS)

$(DATA_DIR):
	mkdir -p $@

$(RESULTS_DIR):
	mkdir -p $@

# Download helper
# $(call download,output,url)
download = curl -L -o $1 $2

$(DATA_DIR)/chr20.1kg.phase3.v5a.vcf.gz: | $(DATA_DIR)
	$(call download,$@,$(URL_chr20))

$(DATA_DIR)/chr21.1kg.phase3.v5a.vcf.gz: | $(DATA_DIR)
	$(call download,$@,$(URL_chr21))

$(DATA_DIR)/example.bed: | $(DATA_DIR)
	$(call download,$@,$(URL_intervals))

datasets: $(DATA_DIR)/chr20.1kg.phase3.v5a.vcf.gz $(DATA_DIR)/chr21.1kg.phase3.v5a.vcf.gz

# Install dependencies (for macOS with Homebrew)
install-deps:
	@echo "Installing benchmark dependencies..."
	@if command -v brew >/dev/null 2>&1; then \
		echo "Installing bcftools and vcftools via Homebrew..."; \
		brew install bcftools vcftools; \
		echo "Installing vcflib via Homebrew..."; \
		brew install vcflib; \
		echo "Installing GATK via Homebrew..."; \
		brew install gatk; \
		echo "All tools installed successfully!"; \
	else \
		echo "Homebrew not found. Please install tools manually:"; \
		echo "  - bcftools: https://samtools.github.io/bcftools/"; \
		echo "  - vcftools: https://vcftools.github.io/"; \
		echo "  - vcflib: https://github.com/vcflib/vcflib"; \
		echo "  - GATK: https://gatk.broadinstitute.org/"; \
	fi

# Tools check
tools-check:
	@echo "Checking VCFX tools..."
	@test -x $(VCFX_VARIANT_COUNTER) || (echo "VCFX_variant_counter not built" && exit 1)
	@test -x $(VCFX_ALLELE_FREQ_CALC) || (echo "VCFX_allele_freq_calc not built" && exit 1)
	@test -x $(VCFX_VARIANT_CLASSIFIER) || (echo "VCFX_variant_classifier not built" && exit 1)
	@test -x $(VCFX_MISSING_DETECTOR) || (echo "VCFX_missing_detector not built" && exit 1)
	@test -x $(VCFX_VALIDATOR) || (echo "VCFX_validator not built" && exit 1)
	@test -x $(VCFX_RECORD_FILTER) || echo "Warning: VCFX_record_filter not built - filtering benchmarks will be skipped"
	@echo "Checking comparison tools..."
	@command -v $(BCFTOOLS) >/dev/null 2>&1 || echo "Warning: bcftools not found - some comparisons will be skipped"
	@command -v $(VCFTOOLS) >/dev/null 2>&1 || echo "Warning: vcftools not found - some comparisons will be skipped"
	@command -v $(GATK) >/dev/null 2>&1 || echo "Warning: GATK not found - some comparisons will be skipped"
	@command -v $(VCFLIB_vcfstats) >/dev/null 2>&1 || echo "Warning: vcflib (vcfstats) not found - some comparisons will be skipped"
	@command -v $(VCFLIB_vcffilter) >/dev/null 2>&1 || echo "Warning: vcflib (vcffilter) not found - some comparisons will be skipped"

# Generate tasks.mk from YAML
tasks.mk: $(MK_TASKS) $(TASKS_FILE)
	$(MK_TASKS) $(TASKS_FILE) > $@

clean:
	rm -rf $(DATA_DIR) $(RESULTS_DIR) tasks.mk $(RESULTS_CSV)
