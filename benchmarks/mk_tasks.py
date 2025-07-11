#!/usr/bin/env python3
"""Generate Makefile rules from tasks.yaml."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml


def main(task_yaml: str) -> None:
    cfg = yaml.safe_load(Path(task_yaml).read_text())
    tasks = cfg.get("tasks", [])

    names = " ".join(t["name"] for t in tasks)
    print(f"TASKS := {names}")
    for t in tasks:
        output = f"results/{t['name']}.csv"
        input_file = f"$(DATA_DIR)/{t['input']}"
        tool = f"$({t['tool']})"
        print()
        print(f"{output}: {input_file}")
        print(f"\t$(call run-task,{tool} {input_file},{output})")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: mk_tasks.py TASKS.yaml")
        raise SystemExit(1)
    main(sys.argv[1])
