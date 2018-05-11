#!/bin/bash
set -x
script_dir="$(dirname $0)"
java -Xms16G -Xmx16G -jar ${script_dir}/trimmomatic.jar "$@"
