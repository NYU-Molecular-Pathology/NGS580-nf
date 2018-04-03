#!/bin/bash
script_dir="$(dirname $0)"
set -x
java -Xms8G -Xmx8G -jar ${script_dir}/GenomeAnalysisTK.jar "$@"
