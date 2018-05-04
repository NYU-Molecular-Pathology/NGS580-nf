Channel.from(
    // [downdb_param, file_name_pattern, protocol_param]
    // ["1000g2015aug", "ALL.sites.2015_08.txt", "1000g2015aug_all"],
    ["refGene", "refGene.txt", "refGene"],
    ["cosmic70", "cosmic70.txt", "cosmic70"]
    // ["clinvar_20170905", "clinvar_20170905.txt", "clinvar_20170905"],
    // ["intervar_20170202", "intervar_20170202.txt", "intervar_20170202"],
    // ["dbnsfp33a", "dbnsfp33a.txt", "dbnsfp33a"]
    // ["esp6500siv2_all", "esp6500siv2_all.txt", "esp6500siv2_all"],
    // ["kaviar_20150923", "kaviar_20150923.txt", "kaviar_20150923"],
    // ["gnomad_exome", "gnomad_exome.txt", "gnomad_exome"],
    // ["gnomad_genome", "gnomad_genome.txt", "gnomad_genome"],
    // ["avsnp150", "avsnp150.txt", "avsnp150"],
    // ["cadd13gt10", "cadd13gt10.txt", ""],
    // ["fathmm", "fathmm.txt", "fathmm"],
    // ["eigen", "eigen.txt", "eigen"]
    )
    .set { annovar_params }
Channel.fromPath("${params.ANNOVAR_DB_DIR}").set { annovar_db_dir }


process make_ANNOVAR_db {
    tag "${downdb_param}-${protocol_param}"
    echo true

    storeDir "${params.ANNOVAR_DB_DIR}"

    input:
    set val(downdb_param), val(output_file), val(protocol_param), file(annovar_db_dir) from annovar_params.combine(annovar_db_dir)

    output:
    file(output_name)

    script:
    output_name = "${params.ANNOVAR_BUILD_VERSION}_${output_file}"
    """
    annotate_variation.pl \
    -downdb \
    -buildver ${params.ANNOVAR_BUILD_VERSION} \
    -webfrom annovar \
    "${downdb_param}" \
    "${annovar_db_dir}"
    """
}
