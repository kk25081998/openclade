#!/bin/sh
# Runs TNT on every script. Usage: TNT=/path/to/tnt sh run_all.sh
cd "$(dirname "$0")"
for f in runs/*.tnt; do echo "$f"; "${TNT:-tnt}" proc "$f"; ; done
