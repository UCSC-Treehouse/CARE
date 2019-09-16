
# Build:
# docker build -f Dockerfile.interactive -t protocol-interactive .
# Run:
# ./interactive.sh dockername port
# (Default port: 18881)
# Once spun up, run these commands in a terminal: (eg for v8):
# ln -s /treehouse/archive/compendium/v8/ /work/cohort
# ln -s /treehouse/archive/downstream/ /work/inputs
# ln -s /treehouse/archive/references/ /work/references
# ln -s /app /work/outputs/app
# ln -s /work/outputs /work/rollup
# Then to run:
# place manifest.tsv in /app
# python /app/run.py

port=${2:-18881}

mkdir outputs
docker run --rm -it -p${port}:8888 \
    -v `pwd`:/app \
    -v `pwd`/outputs:/work/outputs \
    -v /private/groups/treehouse:/treehouse:ro \
    --name $1 protocol-interactive:latest

