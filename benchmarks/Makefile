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

URL_chr20 := https://example.com/chr20_1_1M.vcf.gz
SHA_chr20 := 0000000000000000000000000000000000000000000000000000000000000000
URL_gnomad := https://example.com/gnomad_sv_chr2.vcf.gz
SHA_gnomad := 0000000000000000000000000000000000000000000000000000000000000000
URL_intervals := https://example.com/test_intervals.bed
SHA_intervals := 0000000000000000000000000000000000000000000000000000000000000000

.PHONY: all datasets tools-check

all: datasets tools-check $(RESULTS_CSV)

$(RESULTS_CSV): harness.mk tasks.mk
$(MAKE) -f harness.mk RESULTS_CSV=$(RESULTS_CSV) BENCH_SCRIPTS=$(BENCH_SCRIPTS)

$(DATA_DIR):
mkdir -p $@

$(RESULTS_DIR):
mkdir -p $@

# Download helper
# $(call download,output,url,sha256)
download = wget -O $1 $2 && echo "$3  $1" | sha256sum -c -

$(DATA_DIR)/chr20_1_1M.vcf.gz: | $(DATA_DIR)
$(call download,$@,$(URL_chr20),$(SHA_chr20))

$(DATA_DIR)/gnomad_sv_chr2.vcf.gz: | $(DATA_DIR)
$(call download,$@,$(URL_gnomad),$(SHA_gnomad))

$(DATA_DIR)/test_intervals.bed: | $(DATA_DIR)
$(call download,$@,$(URL_intervals),$(SHA_intervals))

datasets: $(DATA_DIR)/chr20_1_1M.vcf.gz $(DATA_DIR)/gnomad_sv_chr2.vcf.gz $(DATA_DIR)/test_intervals.bed

# Tools check

tools-check:
@test -x $(VCFX_VARIANT_COUNTER) || (echo "VCFX tools not built" && exit 1)

# Generate tasks.mk from YAML

tasks.mk: $(MK_TASKS) $(TASKS_FILE)
$(MK_TASKS) $(TASKS_FILE) > $@

clean:
rm -rf $(DATA_DIR) $(RESULTS_DIR) tasks.mk $(RESULTS_CSV)
