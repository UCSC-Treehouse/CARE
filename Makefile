IMAGE=ucsctreehouse/care:0.15.0.0

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
