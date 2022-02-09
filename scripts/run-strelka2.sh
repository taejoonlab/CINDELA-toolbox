#!/bin/bash
DIR="<directory created by run-conf-strelka2.sh>"

for PY in $(ls $DIR/*py)
do
    python2 $PY -m local
done
