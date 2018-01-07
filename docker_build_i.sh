#!/bin/bash
DOCKER_DIR="docker_build_i"

cp -p tsv2rdf_instruct.py $DOCKER_DIR
cp -p templ_instruct.ttl $DOCKER_DIR
cp -p templ_instruct.ttl.prefix $DOCKER_DIR
cp -p templ_instruct.ttl.evi $DOCKER_DIR
sudo docker build -t med2rdf/instruct $DOCKER_DIR
