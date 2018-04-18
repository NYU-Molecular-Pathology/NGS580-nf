#!/bin/bash
script_dir="$(dirname $0)"
java -jar ${script_dir}/trimmomatic.jar "$@"
# -Xms16G -Xmx16G
