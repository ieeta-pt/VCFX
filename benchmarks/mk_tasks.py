#!/usr/bin/env python3
"""Generate Makefile rules from tasks.yaml."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml


def generate_tool_command(task: dict) -> str:
    """Generate the appropriate command for a task based on the tool type."""
    tool = task["tool"]
    input_file = f"$(DATA_DIR)/{task['input']}"

    if tool == "VCFX_VARIANT_COUNTER":
        return f"$({tool}) < {input_file}"
    elif tool == "VCFX_ALLELE_FREQ_CALC":
        return f"$({tool}) < {input_file}"
    elif tool == "VCFX_VARIANT_CLASSIFIER":
        return f"$({tool}) < {input_file}"
    elif tool == "VCFX_MISSING_DETECTOR":
        return f"$({tool}) < {input_file}"
    elif tool == "VCFX_RECORD_FILTER":
        return f"$({tool}) < {input_file}"
    elif tool == "VCFX_VALIDATOR":
        return f"$(VCFX_VALIDATOR) < {input_file}"
    elif tool == "VCFX_VALIDATOR_MMAP":
        return f"$(VCFX_VALIDATOR) -i {input_file}"
    elif tool == "BCFTOOLS":
        command = task.get("command", "")
        return f"$({tool}) {command} {input_file}"
    elif tool == "VCFTOOLS":
        if "count" in task["name"] or "stress" in task["name"]:
            # vcftools doesn't have --counts, so we'll count sites using --freq
            cmd = f"$({tool}) --vcf {input_file} --freq --stdout 2>/dev/null"
            return f"{cmd} | grep -v '^CHROM' | wc -l"
        elif "allele_freq" in task["name"]:
            return f"$({tool}) --vcf {input_file} --freq --stdout"
        elif "classify" in task["name"]:
            return f"$({tool}) --vcf {input_file} --hardy --stdout"
        elif "missing" in task["name"]:
            return f"$({tool}) --vcf {input_file} --missing-indv --stdout"
        elif "filter" in task["name"]:
            return f"$({tool}) --vcf {input_file} --minQ 30 --recode --stdout"
        else:
            command = task.get("command", "")
            return f"$({tool}) --vcf {input_file} {command} --stdout 2>/dev/null"
    elif tool == "GATK":
        command = task.get("command", "")
        return f"$({tool}) {command} -I {input_file}"
    elif tool == "VCFLIB":
        command = task.get("command", "")
        return f"$({tool}_{command}) {input_file}"
    else:
        # Generic case - assume it's a variable name
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
        elif tool == "VCFLIB":
            # vcflib tools have different names
            vcflib_command = task.get("command", "")
            tool_check = vcflib_command

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
    print("ROOT_DIR := $(abspath ..)")
    print("DATA_DIR := $(ROOT_DIR)/benchmarks/data")
    print("RESULTS_DIR := $(ROOT_DIR)/benchmarks/results")
    print("BENCH_SCRIPTS := $(ROOT_DIR)/benchmarks/scripts")
    print("VCFX_BIN_DIR ?= $(ROOT_DIR)/build/src")
    tools = [
        "VCFX_variant_counter",
        "VCFX_allele_freq_calc",
        "VCFX_variant_classifier",
        "VCFX_missing_detector",
        "VCFX_record_filter",
        "VCFX_validator",
    ]
    for tool in tools:
        var = tool.upper()
        print(f"{var} := $(VCFX_BIN_DIR)/{tool}/{tool}")
    print("BCFTOOLS ?= bcftools")
    print("VCFTOOLS ?= vcftools")
    print("GATK ?= gatk")
    print("VCFLIB_vcfstats ?= vcfstats")
    print("VCFLIB_vcffilter ?= vcffilter")
    print()

    # Generate task list
    names = " ".join(t["name"] for t in tasks)
    print(f"TASKS := {names}")

    # Generate rules for each task
    for task in tasks:
        print(generate_conditional_rule(task))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: mk_tasks.py TASKS.yaml")
        raise SystemExit(1)
    main(sys.argv[1])
