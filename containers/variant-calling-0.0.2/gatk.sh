#!/bin/bash
script_dir="$(dirname $0)"
java -jar ${script_dir}/GenomeAnalysisTK.jar "$@"
# -Xms8G -Xmx8G
