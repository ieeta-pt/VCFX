#!/usr/bin/env python3
"""Generate Makefile rules from tasks.yaml for comprehensive VCFX benchmarking."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml


# Mapping of tool names to their binary paths
VCFX_TOOLS = [
    "VCFX_af_subsetter",
    "VCFX_alignment_checker",
    "VCFX_allele_balance_calc",
    "VCFX_allele_balance_filter",
    "VCFX_allele_counter",
    "VCFX_allele_freq_calc",
    "VCFX_ancestry_assigner",
    "VCFX_ancestry_inferrer",
    "VCFX_annotation_extractor",
    "VCFX_compressor",
    "VCFX_concordance_checker",
    "VCFX_cross_sample_concordance",
    "VCFX_custom_annotator",
    "VCFX_diff_tool",
    "VCFX_distance_calculator",
    "VCFX_dosage_calculator",
    "VCFX_duplicate_remover",
    "VCFX_fasta_converter",
    "VCFX_field_extractor",
    "VCFX_file_splitter",
    "VCFX_format_converter",
    "VCFX_genotype_query",
    "VCFX_gl_filter",
    "VCFX_haplotype_extractor",
    "VCFX_haplotype_phaser",
    "VCFX_header_parser",
    "VCFX_hwe_tester",
    "VCFX_impact_filter",
    "VCFX_inbreeding_calculator",
    "VCFX_indel_normalizer",
    "VCFX_indexer",
    "VCFX_info_aggregator",
    "VCFX_info_parser",
    "VCFX_info_summarizer",
    "VCFX_ld_calculator",
    "VCFX_merger",
    "VCFX_metadata_summarizer",
    "VCFX_missing_data_handler",
    "VCFX_missing_detector",
    "VCFX_multiallelic_splitter",
    "VCFX_nonref_filter",
    "VCFX_outlier_detector",
    "VCFX_phase_checker",
    "VCFX_phase_quality_filter",
    "VCFX_phred_filter",
    "VCFX_population_filter",
    "VCFX_position_subsetter",
    "VCFX_probability_filter",
    "VCFX_quality_adjuster",
    "VCFX_record_filter",
    "VCFX_ref_comparator",
    "VCFX_reformatter",
    "VCFX_region_subsampler",
    "VCFX_sample_extractor",
    "VCFX_sorter",
    "VCFX_subsampler",
    "VCFX_sv_handler",
    "VCFX_validator",
    "VCFX_variant_classifier",
    "VCFX_variant_counter",
]


def generate_tool_command(task: dict) -> str:
    """Generate the appropriate command for a task based on the tool type."""
    tool = task["tool"]
    input_file = f"$(DATA_DIR)/{task['input']}"
    args = task.get("args", "")

    # VCFX tools - most use stdin, validator can use file arg
    if tool == "VCFX_VALIDATOR_MMAP":
        return f"$(VCFX_VALIDATOR) -i {input_file}"
    elif tool.startswith("VCFX_"):
        tool_var = tool.upper()
        if args:
            return f"$({tool_var}) {args} < {input_file}"
        else:
            return f"$({tool_var}) < {input_file}"
    elif tool == "BCFTOOLS":
        command = task.get("command", "")
        return f"$({tool}) {command} {input_file}"
    elif tool == "VCFTOOLS":
        command = task.get("command", "")
        return f"$({tool}) --vcf {input_file} {command} --stdout 2>/dev/null"
    elif tool == "GATK":
        command = task.get("command", "")
        return f"$({tool}) {command} -I {input_file}"
    else:
        # Generic case
        return f"$({tool}) {input_file}"


def generate_conditional_rule(task: dict) -> str:
    """Generate a conditional rule that checks if the tool is available."""
    tool = task["tool"]
    task_name = task["name"]
    input_file = f"$(DATA_DIR)/{task['input']}"
    output_file = f"results/{task_name}.csv"
    command = generate_tool_command(task)

    if tool in ["BCFTOOLS", "VCFTOOLS", "GATK", "VCFLIB"]:
        # For comparison tools, make the rule conditional on tool availability
        tool_check = tool.lower()
        if tool == "GATK":
            tool_check = "gatk"

        return f"""
{output_file}: {input_file}
\t@if command -v {tool_check} >/dev/null 2>&1; then \\
\t\t$(call run-task,{command},{output_file}); \\
\telse \\
\t\techo "Skipping {task_name}: {tool_check} not available" > {output_file}; \\
\t\techo "0" > {output_file}.time; \\
\tfi"""
    else:
        # For VCFX tools, generate normal rule
        return f"""
{output_file}: {input_file}
\t$(call run-task,{command},{output_file})"""


def main(task_yaml: str) -> None:
    cfg = yaml.safe_load(Path(task_yaml).read_text())
    tasks = cfg.get("tasks", [])

    # Generate variable definitions
    print("# Auto-generated from tasks.yaml - do not edit manually")
    print("ROOT_DIR := $(abspath ..)")
    print("DATA_DIR := $(ROOT_DIR)/benchmarks/data")
    print("RESULTS_DIR := $(ROOT_DIR)/benchmarks/results")
    print("BENCH_SCRIPTS := $(ROOT_DIR)/benchmarks/scripts")
    print("VCFX_BIN_DIR ?= $(ROOT_DIR)/build/src")
    print()

    # Generate variable for each VCFX tool
    for tool in VCFX_TOOLS:
        var = tool.upper()
        print(f"{var} := $(VCFX_BIN_DIR)/{tool}/{tool}")

    print()
    print("# Comparison tools (if available)")
    print("BCFTOOLS ?= bcftools")
    print("VCFTOOLS ?= vcftools")
    print("GATK ?= gatk")
    print("VCFLIB_vcfstats ?= vcfstats")
    print("VCFLIB_vcffilter ?= vcffilter")
    print()

    # Generate task list
    names = " ".join(t["name"] for t in tasks)
    print(f"TASKS := {names}")
    print()

    # Generate output file list
    outputs = " ".join(f"results/{t['name']}.csv" for t in tasks)
    print(f"OUTPUTS := {outputs}")
    print()

    # Generate rules for each task
    for task in tasks:
        print(generate_conditional_rule(task))

    print()
    print("# Default target")
    print("all: $(OUTPUTS)")
    print()
    print("# Clean target")
    print("clean:")
    print("\trm -rf results/*.csv results/*.time")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: mk_tasks.py TASKS.yaml")
        raise SystemExit(1)
    main(sys.argv[1])
