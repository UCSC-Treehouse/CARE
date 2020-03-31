IMAGE=ucsctreehouse/care:0.17.1.0
COHORT=v10_polya

run:
	 docker run \
		--rm \
		--user $$UID \
		-v `pwd`/:/work \
		$(IMAGE) run

# adds /app for development
run-dev:
	docker run \
		--rm \
		--user $$UID \
		-v `pwd`/care:/app:ro \
		-v `pwd`/:/work \
		$(IMAGE) run

run-treehouse:
	@echo Running $(IMAGE) vs compendium $(COHORT)
	mkdir -p outputs
	docker run \
		--rm \
		--user $$UID \
		-v `pwd`:/work/rollup:ro \
		-v `pwd`/manifest.tsv:/work/manifest.tsv:ro \
		-v `pwd`/outputs:/work/outputs \
		-v /private/groups/treehouse/archive:/treehouse:ro \
		$(IMAGE) pass-args \
			--inputs /treehouse/downstream \
			--cohort /treehouse/compendium/$(COHORT) \
			--references /treehouse/references
