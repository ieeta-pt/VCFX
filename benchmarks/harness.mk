include tasks.mk

run-task = $(BENCH_SCRIPTS)/run_task.sh "$1" "$2"

results: $(TASKS:%=$(RESULTS_DIR)/%.csv)
	@echo "Results written to $(RESULTS_CSV)"

$(RESULTS_CSV): results
	@echo "task,seconds" > $@
	@for t in $(TASKS); do \
		sec=$$(cat $(RESULTS_DIR)/$$t.csv.time); \
		echo "$$t,$$sec" >> $@; \
	done
