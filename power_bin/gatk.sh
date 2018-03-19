#!/bin/bash
script_dir="$(dirname $0)"
# script_dir="/shared/GATK-P9" # /shared/GATK-P9/GenomeAnalysisTK.jar
export JAVA_TOOL_OPTIONS="-Xmn18g -Xms21g -Xmx21g -XX:ParallelGCThreads=16 -verbose:gc -XX:+PrintGCDetails -XX:+PrintGCTimeStamps -XX:+AggressiveOpts -Djava.library.path=/shared/biobuilds-2017.11/lib/"
set -x
java -Xms16G -Xmx16G -jar ${script_dir}/GenomeAnalysisTK.jar $*
