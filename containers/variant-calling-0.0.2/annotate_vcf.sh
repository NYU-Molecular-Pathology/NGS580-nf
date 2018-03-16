#!/bin/bash

# Annotate a .vcf file with ANNOVAR

# get environment variables
annovar_dir="${ANNOVAR_DIR:-none}" # "/ifs/data/molecpathlab/bin/annovar"
annovar_db_dir="${ANNOVAR_DB_DIR:-none}" # "/ifs/data/molecpathlab/bin/annovar/db/hg19"
annovar_protocol="${ANNOVAR_PROTOCOL:-none}" # "refGene,1000g2015aug_all,clinvar_20170905,intervar_20170202,dbnsfp33a,esp6500siv2_all,kaviar_20150923,gnomad_exome,gnomad_genome,avsnp150,fathmm,eigen"
annovar_operation="${ANNOVAR_OPERATION:-none}" # "g,f,f,f,f,f,f,f,f,f,f,f"
build_version="${ANNOVAR_BUILD_VERSION:-none}" # "hg19"


if [ "${annovar_dir}" != "none" ]; then
    echo "ANNOVAR_DIR is: $annovar_dir"
else
    echo "ERROR: ANNOVAR_DIR is not set, please 'export' it before running"
    exit 1
fi


if [ "${annovar_db_dir}" != "none" ]; then
    echo "ANNOVAR_DB_DIR is: $annovar_db_dir"
else
    echo "ERROR: ANNOVAR_DB_DIR is not set, please 'export' it before running"
    exit 1
fi

if [ "${annovar_protocol}" != "none" ]; then
    echo "ANNOVAR_PROTOCOL is: $annovar_protocol"
else
    echo "ERROR: ANNOVAR_PROTOCOL is not set, please 'export' it before running"
    exit 1
fi

if [ "${annovar_operation}" != "none" ]; then
    echo "ANNOVAR_OPERATION is: $annovar_operation"
else
    echo "ERROR: ANNOVAR_OPERATION is not set, please 'export' it before running"
    exit 1
fi

if [ "${build_version}" != "none" ]; then
    echo "ANNOVAR_BUILD_VERSION is: $build_version"
else
    echo "ERROR: ANNOVAR_BUILD_VERSION is not set, please 'export' it before running"
    exit 1
fi






# get script variable
input_vcf="${1:-none}"
if [ "${input_vcf}" != "none" ]; then
    echo "input_vcf is: $input_vcf"
else
    echo "ERROR: input_vcf was not passed"
    exit 1
fi
[ ! -e "${input_vcf}" ] && echo "ERROR: item does not exist: ${input_vcf}" && exit 1

sample_ID="${2:-fooooooooozzzzz}"
if [ "${sample_ID}" != "fooooooooozzzzz" ]; then
    echo "sample_ID is: $sample_ID"
else
    echo "ERROR: sample_ID was not passed"
    exit 1
fi


# setup output file names
avinput_file="${sample_ID}.avinput"
annovar_output_file="${sample_ID}.${build_version}_multianno.txt"




set -x
# convert to ANNOVAR format
"${annovar_dir}/convert2annovar.pl" --format vcf4old --includeinfo "${input_vcf}" --outfile "${avinput_file}"

# check number of lines between the files
[ ! "$( cat "${avinput_file}" | wc -l )" -eq "$(grep -v '^#' "${input_vcf}" | wc -l)" ] && echo "ERROR: number of entries does not match between files ${input_vcf} and ${avinput_file}" && exit 1 || :

# annovate
"${annovar_dir}/table_annovar.pl" "${avinput_file}" "${annovar_db_dir}" --buildver "${build_version}" --remove --protocol "${annovar_protocol}" --operation "${annovar_operation}" --nastring . --outfile "${sample_ID}"
