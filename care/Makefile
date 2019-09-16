# Set FLAGS in `make run`
# EG make run FLAGS="--clean --stop 2.0"
FLAGS=
PROTOCOL=V13.0.0
COMPENDIUM=v8

build:
	docker build -t protocol:vTEST .

run:
	docker run --user $$UID:`getent group treehouse | cut -d: -f3` --rm -it \
		-v /private/groups/treehouse:/treehouse:ro \
		-v `pwd`/outputs:/work/outputs \
		-v `pwd`/manifest.tsv:/app/manifest.tsv:ro \
		-v `pwd`:/work/rollup:ro \
		protocol:${PROTOCOL} ${FLAGS} \
			--inputs /treehouse/archive/downstream \
			--cohort /treehouse/archive/compendium/${COMPENDIUM} \
			--reference /treehouse/archive/references
			

notebook:
	docker run -it --rm \
		-p 18888:8888 \
		--user root \
    	-e NB_UID=`id -u $$USER` \
    	-e NB_GID=`getent group treehouse | cut -d: -f3` \
		-v `pwd`/:/app \
		-v `pwd`/outputs/:/outputs \
		-v `pwd`:/rollup:ro \
		-v /private/groups/treehouse:/treehouse:ro \
		protocol:${PROTOCOL} -j
