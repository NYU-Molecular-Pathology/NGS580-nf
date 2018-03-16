#!/bin/bash
script_dir="$(dirname $0)"
set -x
java -Xms16G -Xmx16G -jar ${script_dir}/GenomeAnalysisTK.jar $*
