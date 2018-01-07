#!/bin/bash
sudo docker run -v `pwd`:/mnt med2rdf/instruct -c /mnt/tsv2rdf_instruct.json
