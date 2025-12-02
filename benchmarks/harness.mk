include tasks.mk

run-task = $(BENCH_SCRIPTS)/run_task_detailed.sh "$1" "$2"
RESULTS_CSV := results.csv

results: $(TASKS:%=$(RESULTS_DIR)/%.csv)
	@echo "Results written to $(RESULTS_CSV)"

$(RESULTS_CSV): results
	@echo "Generating benchmark summary..."
	@echo "task,tool,input,seconds,status" > $@
	@for t in $(TASKS); do \
		if [ -f $(RESULTS_DIR)/$$t.csv.time ]; then \
			sec=$$(cat $(RESULTS_DIR)/$$t.csv.time); \
			tool=$$(echo $$t | cut -d'_' -f1); \
			input=$$(echo $$t | cut -d'_' -f2-); \
			if [ -s $(RESULTS_DIR)/$$t.csv ] && [ "$$sec" != "0" ]; then \
				echo "$$t,$$tool,$$input,$$sec,completed" >> $@; \
			else \
				echo "$$t,$$tool,$$input,$$sec,skipped" >> $@; \
			fi; \
		else \
			echo "$$t,unknown,unknown,0,failed" >> $@; \
		fi; \
	done
	@echo "Benchmark summary:"
	@echo "=================="
	@cat $@
	@echo ""
	@echo "Performance comparison:"
	@echo "======================="
	@$(BENCH_SCRIPTS)/compare_results.py $@ || echo "Comparison script not available"

comparison: $(RESULTS_CSV)
	@echo "Generating detailed comparison report..."
	@$(BENCH_SCRIPTS)/compare_results.py $(RESULTS_CSV) --detailed || echo "Comparison script not available"
