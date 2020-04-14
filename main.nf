// NGS580 Target exome analysis for 580 gene panel
import java.nio.file.Files;
import java.text.SimpleDateFormat;
import groovy.json.JsonSlurper;
def jsonSlurper = new JsonSlurper()
def workflowTimestamp = "${workflow.start.format('yyyy-MM-dd-HH-mm-ss')}"
def username = System.getProperty("user.name")
String localhostname = java.net.InetAddress.getLocalHost().getHostName();
Date now = new Date()
SimpleDateFormat timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss")
def currentDirPath = new File(System.getProperty("user.dir")).getCanonicalPath()
Process uname_proc = 'uname'.execute()
def uname = "${uname_proc.text.trim()}"
def NXF_PIDFILE = ".nextflow.pid"
if(uname == "Linux"){
    // /proc/self only works on Linux
    int pid = Integer.parseInt(new File("/proc/self").getCanonicalFile().getName())
    File pid_file = new File("${NXF_PIDFILE}")
    pid_file.write("${pid}\n")

}

// ~~~~~~~~~~ CONFIGURATION ~~~~~~~~~~ //
// configure pipeline settings
// order of evaluation
// 1. uses CLI passed arg in 'params'
// 2. uses nextflow.config arg in 'params'
// 3. uses values from config.json
// 4. uses default value

// default values; overriden by nextflow.config and CLI
params.configFile = "config.json"
params.reportDir = "report"
params.outputDir = "output"
params.ref_dir = "ref"
params.nxf_completion_log = ".nextflow.completion.log"
params.GIT_CURRENT_BRANCH = "none"
params.GIT_CURRENT_COMMIT = "none"
params.GIT_CURRENT_TAG = "none"
params.GIT_RECENT_TAG = "none"

// default config values
def defaultParams = [:]
defaultParams.runID = "NGS580_run"
defaultParams.workflowLabel = "NGS580"
defaultParams.numTargetSplitLines = 50
defaultParams.targetsBed = "targets/targets.580.bed"
defaultParams.targetsAnnotatedBed = "targets/targets.annotated.580.bed"
defaultParams.projectDir = "${workflow.projectDir}"
defaultParams.samplesheet = "samples.analysis.tsv"
defaultParams.SeraCareSelectedTsv = "data/SeraCare-selected-variants.tsv"
defaultParams.SeraCareErrorRate = 0.02
defaultParams.CNVPool = "ref/CNV-Pool/CNV-Pool.580.cnn"
defaultParams.HapMapBam = "ref/HapMap-Pool/HapMap-pool.bam"
defaultParams.HapMapBai = "ref/HapMap-Pool/HapMap-pool.bam.bai"

// load the JSON config, if present
def externalConfig
def externalConfigFile_obj = new File("${params.configFile}")
if ( externalConfigFile_obj.exists() ) {
    log.info("Loading configs from ${params.configFile}")
    String externalConfigJSON = externalConfigFile_obj.text
    externalConfig = jsonSlurper.parseText(externalConfigJSON)
}
// add and overwrite all values in defaultParams with items from config
if( externalConfig ){
    log.info("Updating with configs from ${params.configFile}")
    // overwrite existing default params with JSON config values
    defaultParams.keySet().each { key ->
        if( externalConfig.containsKey(key) ){
            defaultParams[key] = externalConfig[key]
        }
    }
    // append the missing JSON configs to default params; doesn't overwrite
    externalConfig.keySet().each { key ->
        if( ! defaultParams.containsKey(key) ){
            defaultParams[key] = externalConfig[key]
        }
    }
}

// add all the entries to params; wont overwrite existing entries from CLI or nextflow.config
params << defaultParams

// set variables to use throughout the pipeline
def numTargetSplitLines = params.numTargetSplitLines
def HapMapBam = params.HapMapBam
def HapMapBai = params.HapMapBai
def targetsBed = params.targetsBed
def targetsAnnotatedBed = params.targetsAnnotatedBed
def runID = params.runID
def samplesheet = params.samplesheet
def SeraCareSelectedTsv = params.SeraCareSelectedTsv
def SeraCareSelectedTsvFile = new File("${SeraCareSelectedTsv}").getName()
def SeraCareErrorRate = params.SeraCareErrorRate
def CNVPool = params.CNVPool
def projectDir = params.runID

// Enable or disable some pipeline steps here TODO: better config management for this
disable_multiqc = true // for faster testing of the rest of the pipeline
disable_msisensor = true // breaks on very small demo datasets
disable_delly2 = true
disable_eval_pair_vcf = true
disable_pindel = true
//disable_varscan2 = true

// load a mapping dict to use for keeping track of the names and suffixes for some files throughout the pipeline
String filemapJSON = new File("filemap.json").text
def filemap = jsonSlurper.parseText(filemapJSON)

// names of some important output files to use throughout the pipeline
def all_annotations_file = filemap.files.all_annotations_file
def all_seracare_annotations_file = filemap.files.all_seracare_annotations_file
def all_HapMapPool_annotations_file = filemap.files.all_HapMapPool_annotations_file
def samplesheet_output_file = filemap.files.samplesheet_output_file
def sample_coverage_file = filemap.files.sample_coverage_file
def interval_coverage_file = filemap.files.interval_coverage_file
def signatures_weights_file = filemap.files.signatures_weights_file
def targets_annotations_file = filemap.files.targets_annotations_file
def tmb_file = filemap.files.tmb_file
def git_json = filemap.files.git_json
def snp_overlap_file = filemap.files.snp_overlap_file
def outputDirPath = new File(params.outputDir).getCanonicalPath()
def reportDirPath = new File(params.reportDir).getCanonicalPath()

// REFERENCE FILES
params.ref_fa = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.ref_fa_bwa_dir = "${params.ref_dir}/BWA/hg19"
params.ref_fai = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
params.ref_dict = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
params.ref_chrom_sizes = "${params.ref_dir}/Illumina/hg19/chrom.sizes"
params.microsatellites = "${params.ref_dir}/msisensor/hg19/microsatellites.list"
params.trimmomatic_contaminant_fa = "${params.ref_dir}/contaminants/trimmomatic.fa"
params.gnomAD_sites_vcf = "${params.gnomAD_dir}/gnomad.exomes.r2.0.2.sites.vcf.gz"
params.gnomAD_sites_tbi = "${params.gnomAD_dir}/gnomad.exomes.r2.0.2.sites.vcf.gz.tbi"
params.vep_cache_dir = "/gpfs/data/molecpathlab/ref/vep"
params.ExAC_sites_vcf = "${params.ExAC_dir}/ExAC.r0.3.sites.vep.hg19.vcf.gz"
params.ExAC_sites_tbi = "${params.ExAC_dir}/ExAC.r0.3.sites.vep.hg19.vcf.gz.tbi"
params.ExAC_dir = "/gpfs/data/molecpathlab/ref/ExAC"
params.gatk_bundle_dir = "${params.ref_dir}/gatk-bundle"
params.gatk_1000G_phase1_indels_hg19_vcf = "${params.gatk_bundle_dir}/1000G_phase1.indels.hg19.vcf"
params.gatk_1000G_phase1_indels_hg19_vcf_idx = "${params.gatk_bundle_dir}/1000G_phase1.indels.hg19.vcf.idx"
params.mills_and_1000G_gold_standard_indels_hg19_vcf = "${params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
params.mills_and_1000G_gold_standard_indels_hg19_vcf_idx = "${params.gatk_bundle_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf.idx"
params.dbsnp_ref_vcf = "${params.gatk_bundle_dir}/dbsnp_138.hg19.vcf"
params.dbsnp_ref_vcf_idx = "${params.gatk_bundle_dir}/dbsnp_138.hg19.vcf.idx"
params.dbsnp_ref_vcf_gz = "${params.gatk_bundle_dir}/dbsnp_138.hg19.vcf.gz"
params.dbsnp_ref_vcf_gz_tbi = "${params.gatk_bundle_dir}/dbsnp_138.hg19.vcf.gz.tbi"
params.cosmic_ref_vcf = "${params.ref_dir}/hg19/CosmicCodingMuts_v73.hg19.vcf"
params.cosmic_ref_vcf_idx = "${params.ref_dir}/hg19/CosmicCodingMuts_v73.hg19.vcf.idx"
params.common_snp_vcf = "${params.ref_dir}/hg19/common_all_20170710.vcf.gz"
params.common_snp_vcf_tbi = "${params.ref_dir}/hg19/common_all_20170710.vcf.gz.tbi"
params.germline_resource_gz = "/gpfs/data/molecpathlab/ref/gatk-bundle/af-only-gnomad.raw.sites.hg19.vcf.gz" // added for gatk4
params.germline_resource_gz_tbi = "/gpfs/data/molecpathlab/ref/gatk-bundle/af-only-gnomad.raw.sites.hg19.vcf.gz.tbi" // added for gatk4

// ~~~~~~~~~~ VALIDATION ~~~~~~~~~~ //
// make sure reference data directory exists
def ref_dir = new File("${params.ref_dir}")
if( !ref_dir.exists() ){
    log.error "Ref dir does not exist: ${params.ref_dir}"
    exit 1
}

// make sure annovar reference databases directory exist
def ANNOVAR_DB_DIR = new File("${params.ANNOVAR_DB_DIR}")
if( !ANNOVAR_DB_DIR.exists() ){
    log.error "ANNOVAR database dir does not exist: ${params.ANNOVAR_DB_DIR}"
    exit 1
}

// ~~~~~ START WORKFLOW ~~~~~ //
log.info "~~~~~~~ NGS580 Pipeline ~~~~~~~"
if(uname == "Linux"){
    // /proc/self only works on Linux
    int pid = Integer.parseInt(new File("/proc/self").getCanonicalFile().getName())
    log.info "* pid:                ${pid}"
}
log.info "* hostname:           ${localhostname}"
log.info "* uname:              ${uname}"
log.info "* Launch time:        ${workflowTimestamp}"
log.info "* Run ID:             ${runID}"
log.info "* Samplesheet:        ${samplesheet}"
log.info "* Targets:            ${targetsBed}"
log.info "* Project dir:        ${workflow.projectDir}"
log.info "* Launch dir:         ${workflow.launchDir}"
log.info "* Work dir:           ${workflow.workDir.toUriString()}"
log.info "* Output dir:         ${outputDirPath}"
log.info "* Profile:            ${workflow.profile ?: '-'}"
log.info "* Script name:        ${workflow.scriptName ?: '-'}"
log.info "* Script ID:          ${workflow.scriptId ?: '-'}"
log.info "* Container engine:   ${workflow.containerEngine?:'-'}"
log.info "* Workflow session:   ${workflow.sessionId}"
log.info "* Nextflow run name:  ${workflow.runName}"
log.info "* Nextflow version:   ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})"
log.info "* Launch command:\n${workflow.commandLine}\n"

// ~~~~~ DATA INPUT ~~~~~ //
// targets .bed file
Channel.fromPath( file(targetsBed) ).set{ targets_bed } // TODO: why is this here? duplicated..
Channel.fromPath( file(targetsAnnotatedBed) ).set{ targets_annotated_bed }

// reference files
Channel.fromPath( file(targetsBed) ).into { targets_bed;
    targets_bed2;
    targets_bed3;
    targets_bed4;
    targets_bed5;
    targets_bed6;
    targets_bed7;
    targets_bed8;
    targets_bed9;
    targets_bed10;
    targets_bed11;
    targets_bed12;
    targets_bed13;
    targets_bed14;
    targets_bed15
  }

Channel.fromPath( file(params.ref_fa) ).into { ref_fasta;
    ref_fasta2;
    ref_fasta3;
    ref_fasta4;
    ref_fasta5;
    ref_fasta6;
    ref_fasta7;
    ref_fasta8;
    ref_fasta9;
    ref_fasta10;
    ref_fasta11;
    ref_fasta12;
    ref_fasta13;
    ref_fasta14;
    ref_fasta15;
    ref_fasta16;
    ref_fasta17;
    ref_fasta18;
    ref_fasta19;
    ref_fasta20;
    ref_fasta21;
    ref_fasta22;
    ref_fasta23;
    ref_fasta24
  }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai;
    ref_fai2;
    ref_fai3;
    ref_fai4;
    ref_fai5;
    ref_fai6;
    ref_fai7;
    ref_fai8;
    ref_fai9;
    ref_fai10;
    ref_fai11;
    ref_fai12;
    ref_fai13;
    ref_fai14;
    ref_fai15;
    ref_fai16;
    ref_fai17;
    ref_fai18;
    ref_fai19;
    ref_fai20;
    ref_fai21;
    ref_fai22;
    ref_fai23;
    ref_fai24
  }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict;
    ref_dict2;
    ref_dict3;
    ref_dict4;
    ref_dict5;
    ref_dict6;
    ref_dict7;
    ref_dict8;
    ref_dict9;
    ref_dict10;
    ref_dict11;
    ref_dict12;
    ref_dict13;
    ref_dict14;
    ref_dict15;
    ref_dict16;
    ref_dict17;
    ref_dict18;
    ref_dict19;
    ref_dict20;
    ref_dict21;
    ref_dict22
  }

Channel.fromPath( file(params.ref_chrom_sizes) ).set{ ref_chrom_sizes }
Channel.fromPath( file(params.trimmomatic_contaminant_fa) ).set{ trimmomatic_contaminant_fa }
Channel.fromPath( file(params.ref_fa_bwa_dir) ).set{ ref_fa_bwa_dir }

Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf) ).into{ gatk_1000G_phase1_indels_vcf;
    gatk_1000G_phase1_indels_vcf2;
    gatk_1000G_phase1_indels_vcf3;
    gatk_1000G_phase1_indels_vcf4 }
Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf_idx) ).into{ gatk_1000G_phase1_indels_vcf_idx;
    gatk_1000G_phase1_indels_vcf_idx2;
    gatk_1000G_phase1_indels_vcf_idx3;
    gatk_1000G_phase1_indels_vcf_idx4 }

Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf) ).into{ mills_and_1000G_gold_standard_indels_vcf;
    mills_and_1000G_gold_standard_indels_vcf2;
    mills_and_1000G_gold_standard_indels_vcf3;
    mills_and_1000G_gold_standard_indels_vcf4 }
Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf_idx) ).into{ mills_and_1000G_gold_standard_indels_vcf_idx;
    mills_and_1000G_gold_standard_indels_vcf_idx2;
    mills_and_1000G_gold_standard_indels_vcf_idx3;
    mills_and_1000G_gold_standard_indels_vcf_idx4 }

Channel.fromPath( file(params.dbsnp_ref_vcf) ).into{ dbsnp_ref_vcf;
    dbsnp_ref_vcf2;
    dbsnp_ref_vcf3;
    dbsnp_ref_vcf4;
    dbsnp_ref_vcf5;
    dbsnp_ref_vcf6;
    dbsnp_ref_vcf7;
    dbsnp_ref_vcf8 }
Channel.fromPath( file(params.dbsnp_ref_vcf_idx) ).into{ dbsnp_ref_vcf_idx;
    dbsnp_ref_vcf_idx2;
    dbsnp_ref_vcf_idx3;
    dbsnp_ref_vcf_idx4;
    dbsnp_ref_vcf_idx5;
    dbsnp_ref_vcf_idx6;
    dbsnp_ref_vcf_idx7;
    dbsnp_ref_vcf_idx8 }
Channel.fromPath( file(params.dbsnp_ref_vcf_gz) ).set { dbsnp_ref_vcf_gz }
Channel.fromPath( file(params.dbsnp_ref_vcf_gz_tbi) ).set { dbsnp_ref_vcf_gz_tbi }

Channel.fromPath( file(params.germline_resource_gz) ).set { germline_resource_gz } //added for gatk4
Channel.fromPath( file(params.germline_resource_gz_tbi) ).set { germline_resource_gz_tbi }//added for gatk4

Channel.fromPath( file(params.cosmic_ref_vcf) ).into{ cosmic_ref_vcf; cosmic_ref_vcf2 }
Channel.fromPath( file(params.cosmic_ref_vcf_idx) ).into{ cosmic_ref_vcf_idx; cosmic_ref_vcf_idx2 }

Channel.fromPath( file(params.common_snp_vcf) ).set{ common_snp_vcf }
Channel.fromPath( file(params.common_snp_vcf_tbi) ).set{ common_snp_vcf_tbi }

Channel.fromPath( file(params.microsatellites) ).set{ microsatellites }
Channel.fromPath( file(params.ANNOVAR_DB_DIR) ).into { annovar_db_dir;
    annovar_db_dir2;
    annovar_db_dir3;
    annovar_db_dir4 }

// report and output dir
Channel.fromPath("${outputDirPath}").into { analysis_output; analysis_output2 }

// load analysis report files
Channel.fromPath("${reportDirPath}/util").into { report_utils; report_utils2 }

Channel.fromPath("${reportDirPath}/analysis/*")
        .set { analysis_report_files_base }
analysis_report_files_base.mix(report_utils).set { analysis_report_files }

// load samples report files
Channel.fromPath("${reportDirPath}/samples/*")
        .toList()
        .map { items ->
            return( [items])
        }
        .set { samples_report_files }

// read samples from analysis samplesheet
Channel.fromPath( file(samplesheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sampleID = row['Sample']
            def reads1 = row['R1'].tokenize( ',' ).collect { file(it) } // comma-sep string into list of files
            def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
            return [ sampleID, reads1, reads2 ]
        }
        .tap { samples_R1_R2; samples_R1_R2_2 } // set of all fastq R1 R2 per sample
        .map { sampleID, reads1, reads2 ->
            return [ reads1, reads2 ]
        }
        .flatMap().flatMap()
        .set { samples_each_fastq } // emit each fastq file individually, no sampleID


// read sample tumor-normal pairs from analysis sheet
Channel.fromPath( file(samplesheet) )
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def tumorID = row['Tumor']
            def normalID = row['Normal']
            return [ tumorID, normalID ]
        }
        .filter { item ->
            item[0] != item[1] // remove Normal samples
        }
        .tap { samples_pairs_with_NA; samples_pairs_with_NA2 }
        .filter { item ->
            item[1] != 'NA' // unpaired samples
        }
        .into { samples_pairs;
        samples_pairs2;
        samples_pairs3;
        samples_pairs4 }

// read sample IDs from analysis sheet
Channel.fromPath( file(samplesheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sampleID = row['Sample']
            return(sampleID)
        }
        .unique()
        .into { sampleIDs; sampleIDs2; sampleIDs3 }

// find all the HapMap samples
sampleIDs2.filter {
    def is_hapmap = "${it}".toLowerCase().contains("hapmap")
    return(is_hapmap)
}.set { hapmap_sample_ids }

// ref HapMap Pool bam and bai
Channel.from([ [file("${HapMapBam}"), file("${HapMapBai}")] ]).filter { items ->
    def bam = items[0]
    def bai = items[1]
    // only run these steps if params were set and files exist
    HapMapBam && HapMapBai && bam.exists() && bai.exists()
    }.set { hapmap_pool_ch }

// find all the SeraCare samples
sampleIDs3.filter {
    def is_seracare = "${it}".toLowerCase().contains("seracare")
    return(is_seracare)
}.set { seracare_sample_ids }

Channel.fromPath( "${SeraCareSelectedTsv}").into { seracare_selected_tsv; seracare_selected_tsv2 }
Channel.fromPath("${CNVPool}").set { cnv_pool_ch }
Channel.fromPath( file(samplesheet) ).set { samples_analysis_sheet }

// reference file for VEP calling
Channel.fromPath( file(params.gnomAD_sites_vcf) ).set{ gnomAD_sites_vcf}
Channel.fromPath( file(params.gnomAD_sites_tbi) ).set{ gnomAD_sites_tbi}
Channel.fromPath( file(params.ExAC_sites_vcf) ).set{ ExAC_sites_vcf}
Channel.fromPath( file(params.ExAC_sites_tbi) ).set{ ExAC_sites_tbi}
Channel.fromPath( file(params.vep_cache_dir) ).set{ vep_cache_dir }

// logging channels
Channel.from("Sample\tProgram\tType\tNote\tFiles").set { failed_samples }
Channel.from("Comparison\tTumor\tNormal\tChunk\tProgram\tProgramType\tNote\tFiles").set { failed_pairs }

// ~~~~~ PIPELINE TASKS ~~~~~ //
// PREPROCESSING
process git {
    // get repo information
    publishDir "${params.outputDir}", mode: 'copy'
    output:
    file("${output_file}") into git_json_ch

    script:
    output_file = "${git_json}"
    """
    git.py --dir "${workflow.projectDir}" -o "${output_file}"
    """
}

process copy_samplesheet {
    // make a copy of the samplesheet in the output directory
    // this ensures the output sheet has the correct name
    publishDir "${params.outputDir}", mode: 'copy'
    executor "local"

    input:
    file(input_sheet: "input_samplesheet.tsv") from samples_analysis_sheet

    output:
    file("${samplesheet_output_file}") into samplesheet_output_file_ch
    val("") into done_copy_samplesheet

    script:
    """
    cp "${input_sheet}" "${samplesheet_output_file}"
    """
}


process print_metadata {
    // print the workflow meta data to the output directory
    publishDir "${params.outputDir}", mode: 'copy'
    executor "local"

    input:
    val(x) from Channel.from('')

    output:
    file("${output_file}") into metadata_ch
    val("") into done_print_metadata

    script:
    output_file = "meta.tsv"
    """
    printf "Run\tTime\tSession\tWorkflow\tLocation\tSystem\tOutputPath\tUsername\n" > "${output_file}"
    printf "${runID}\t${workflowTimestamp}\t${workflow.sessionId}\t${workflow.runName}\t${workflow.projectDir}\t${localhostname}\t${outputDirPath}\t${username}\n" >> "${output_file}"
    """
}

process targets_zip {
    input:
    file(targets_bed) from targets_bed13

    output:
    set file("${output_bgz}"), file("${output_index}") into targets_zipped, targets_zipped2

    script:
    output_bgz = "targets.bed.bgz"
    output_index = "targets.bed.bgz.tbi"
    """
    sort -V -k1,1 -k2,2 "${targets_bed}" > targets.sorted.bed
    bgzip -c targets.sorted.bed > "${output_bgz}"
    tabix -p bed "${output_bgz}"
    """
}

process targets_metrics {
    // print metrics about the targets
    publishDir "${params.outputDir}/metrics/targets", mode: 'copy', pattern: "*.bed"
    publishDir "${params.outputDir}", mode: 'copy', pattern: "*${targets_metrics}"
    executor "local"

    input:
    file(targets) from targets_bed9

    output:
    file("${targets_metrics}") into targets_metrics_ch

    script:
    targets_sorted = "targets.sorted.bed"
    targets_merged = "targets.merged.bed"
    targets_metrics = "targets.metrics.tsv"
    """
    num_targets="\$(cat "${targets}" | wc -l)"
    targets_md5="\$(python -c "import hashlib; print(hashlib.md5(open('${targets}', 'rb').read()).hexdigest())")"

    # check if there are strands in the targets
    if [ "\$(bed.py "${targets}" hasStrands)" == "True" ]; then
        sort -k 1,1 -k2,2n "${targets}" > "${targets_sorted}"
        bedtools merge -s -i "${targets_sorted}" > "${targets_merged}"
    else
        sort -k 1,1 -k2,2n "${targets}" > "${targets_sorted}"
        bedtools merge -i "${targets_sorted}" > "${targets_merged}"
    fi

    num_merged_targets="\$(cat "${targets_merged}" | wc -l)"

    targets_coverage_bp="\$(bed.py "${targets}" breadthOfCoverage)"
    targets_coverage_Mbp="\$(python -c "print( \${targets_coverage_bp} / float((10**6)) )")"

    targets_filename="\$(python -c "import os; print(os.path.basename(os.path.realpath('${targets}')))")"

    printf 'Targets File\tNumber of Targets\tNumber of Merged Targets\tBreadth Of Coverage (Mbp)\tBreadth Of Coverage (bp)\tmd5\n' > "${targets_metrics}"
    printf "\${targets_filename}\t\${num_targets}\t\${num_merged_targets}\t\${targets_coverage_Mbp}\t\${targets_coverage_bp}\t\${targets_md5}\n" >> "${targets_metrics}"
    """
}

target_ANNOVAR_BUILD_VERSION = "hg19"
target_ANNOVAR_PROTOCOL = "refGene"
target_ANNOVAR_OPERATION = "g"
process annotate_targets {
    // annotate the target.bed regions
    input:
    set file(targets_bed), file(annovar_db_dir) from targets_bed11.combine(annovar_db_dir4)

    output:
    file("${output_file}") into annotated_targets

    script:
    prefix = "targets"
    interval_tmp = "${prefix}.intervals.tmp"
    remainder_tsv = "${prefix}.remainder.tsv"
    avinput_file = "${prefix}.avinput"
    annovar_output_txt = "${prefix}.${target_ANNOVAR_BUILD_VERSION}_multianno.txt"
    output_file = "${targets_annotations_file}"
    """
    # convert table to ANNOVAR format for annotation; http://annovar.openbioinformatics.org/en/latest/user-guide/input/
    # add '0' cols for ref and alt
    cut -f1-3 "${targets_bed}" > "${avinput_file}"
    sed -e 's|\$|\t0|' -i "${avinput_file}"
    sed -e 's|\$|\t0|' -i "${avinput_file}"

    # annotate
    table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
    --buildver "${target_ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${target_ANNOVAR_PROTOCOL}" \
    --operation "${target_ANNOVAR_OPERATION}" \
    --nastring . \
    --outfile "${prefix}"

    mv "${annovar_output_txt}" "${output_file}"
    """
}

process fastq_merge {
    // merge multiple R1 and R2 fastq files (e.g. split by lane) into a single fastq each
    input:
    set val(sampleID), file(fastq_r1: "*"), file(fastq_r2: "*") from samples_R1_R2

    output:
    set val(sampleID), file("${merged_fastq_R1}"), file("${merged_fastq_R2}") into samples_fastq_merged
    file("${num_reads_R1}")
    file("${num_reads_R2}")
    file("${num_reads}")
    val(sampleID) into done_fastq_merge

    script:
    prefix = "${sampleID}"
    merged_fastq_R1 = "${prefix}_R1.fastq.gz"
    merged_fastq_R2 = "${prefix}_R2.fastq.gz"
    num_reads = "${prefix}.reads.txt"
    num_reads_R1 = "${prefix}_R1.reads.txt"
    num_reads_R2 = "${prefix}_R2.reads.txt"
    """
    cat ${fastq_r1} > "${merged_fastq_R1}"
    cat ${fastq_r2} > "${merged_fastq_R2}"

    # get the number of reads
    zcat "${merged_fastq_R1}" | awk '{s++}END{print s/4}' > "${num_reads}"
    cp "${num_reads}" "${num_reads_R1}"
    zcat "${merged_fastq_R2}" | awk '{s++}END{print s/4}' > "${num_reads_R2}"
    """
}


process trimmomatic {
    // trim low quality bases from reads
    publishDir "${params.outputDir}/reads/trimmed", mode: 'copy', pattern: "*${fastq_R1_trimmed}"
    publishDir "${params.outputDir}/reads/trimmed", mode: 'copy', pattern: "*${fastq_R2_trimmed}"
    publishDir "${params.outputDir}/reads/stats", mode: 'copy', pattern: "*.txt"

    input:
    set val(sampleID), file(read1), file(read2), file(trimmomatic_contaminant_fa) from samples_fastq_merged.combine(trimmomatic_contaminant_fa)

    output:
    set val(sampleID), file("${fastq_R1_trimmed}"), file("${fastq_R2_trimmed}") into samples_fastq_trimmed, samples_fastq_trimmed2
    file("${num_reads_trim}")
    file("${num_reads_trim_R1}")
    file("${num_reads_trim_R2}")
    file("${num_reads_unpaired_R1}")
    file("${num_reads_unpaired_R2}")
    val(sampleID) into done_trimmomatic

    script:
    prefix = "${sampleID}"
    fastq_R1_trimmed = "${prefix}_R1.trim.fastq.gz"
    fastq_R2_trimmed = "${prefix}_R2.trim.fastq.gz"
    fastq_R1_unpaired = "${prefix}_R1.unpaired.fastq.gz"
    fastq_R2_unpaired = "${prefix}_R2.unpaired.fastq.gz"
    num_reads_trim = "${prefix}.trim.reads.txt"
    num_reads_trim_R1 = "${prefix}_R1.trim.reads.txt"
    num_reads_trim_R2 = "${prefix}_R2.trim.reads.txt"
    num_reads_unpaired_R1 = "${prefix}_R1.unpaired.reads.txt"
    num_reads_unpaired_R2 = "${prefix}_R2.unpaired.reads.txt"
    """
    trimmomatic.sh PE -threads \${NSLOTS:-\${NTHREADS:-1}} \
    "${read1}" "${read2}" \
    "${fastq_R1_trimmed}" "${fastq_R1_unpaired}" \
    "${fastq_R2_trimmed}" "${fastq_R2_unpaired}" \
    ILLUMINACLIP:${trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35

    # get the number of reads
    zcat "${fastq_R1_trimmed}" | awk '{s++}END{print s/4}' > "${num_reads_trim}"
    cp "${num_reads_trim}" "${num_reads_trim_R1}"
    zcat "${fastq_R2_trimmed}" | awk '{s++}END{print s/4}' > "${num_reads_trim_R2}"
    zcat "${fastq_R1_unpaired}" | awk '{s++}END{print s/4}' > "${num_reads_unpaired_R1}"
    zcat "${fastq_R2_unpaired}" | awk '{s++}END{print s/4}' > "${num_reads_unpaired_R2}"
    """
}


process fastqc {
    // quality control checking with FastQC
    publishDir "${params.outputDir}/qc/fastqc", mode: 'copy'

    input:
    set val(sampleID),  file(fastq_R1), file(fastq_R2) from samples_fastq_trimmed2

    output:
    file(output_R1_html)
    file(output_R1_zip)
    file(output_R2_html)
    file(output_R2_zip)
    val(sampleID) into done_fastqc_trim

    script:
    output_R1_html = "${fastq_R1}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
    output_R1_zip = "${fastq_R1}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
    output_R2_html = "${fastq_R2}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
    output_R2_zip = "${fastq_R2}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
    """
    fastqc -o . "${fastq_R1}"
    fastqc -o . "${fastq_R2}"
    """

}

process alignment {
    // first pass alignment with BWA
    publishDir "${params.outputDir}/alignments/raw", mode: 'copy'

    input:
    set val(sampleID), file(fastq_R1_trim), file(fastq_R2_trim), file(ref_fa_bwa_dir) from samples_fastq_trimmed.combine(ref_fa_bwa_dir)

    output:
    set val(sampleID), file("${bam_file}") into samples_bam, samples_bam2, samples_bam3, samples_bam4
    val(sampleID) into done_alignment

    script:
    prefix = "${sampleID}"
    bam_file = "${prefix}.bam"
    """
    bwa mem \
    -M -v 1 \
    -t \${NSLOTS:-\${NTHREADS:-1}} \
    -R '@RG\\tID:${sampleID}\\tSM:${sampleID}\\tLB:${sampleID}\\tPL:ILLUMINA' \
    "${ref_fa_bwa_dir}/genome.fa" \
    "${fastq_R1_trim}" "${fastq_R2_trim}" | \
    sambamba view \
    --sam-input \
    --nthreads=\${NSLOTS:-\${NTHREADS:-1}} \
    --filter='mapping_quality>=10' \
    --format=bam \
    --compression-level=0 \
    /dev/stdin | \
    sambamba sort \
    --nthreads=\${NSLOTS:-\${NTHREADS:-1}} \
    --out="${bam_file}" /dev/stdin
    """
}


process samtools_flagstat {
    // alignment stats
    publishDir "${params.outputDir}/alignments/stats", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam) from samples_bam

    output:
    set val(sampleID), file("${flagstat}") into flagstats
    val(sampleID) into done_sambamba_flagstat

    script:
    prefix = "${sampleID}"
    flagstat = "${prefix}.flagstat.txt"
    """
    samtools flagstat "${sample_bam}" > "${flagstat}"
    """
}

process samtools_flagstat_table {
    // convert flagstat output to a flat table
    publishDir "${params.outputDir}/alignments/stats", mode: 'copy'

    input:
    set val(sampleID), file(flagstat) from flagstats

    output:
    file("${output_file}") into sambamba_flagstat_tables
    val(sampleID) into done_sambamba_flagstat_table

    script:
    prefix = "${sampleID}"
    output_file = "${prefix}.flagstat.tsv"
    """
    flagstat2table.R "${flagstat}" tmp.tsv
    paste-col.py -i tmp.tsv --header "Sample" -v "${sampleID}" > "${output_file}"
    """
}
sambamba_flagstat_tables.collectFile(name: ".flagstat.tsv", keepHeader: true).set { sambamba_flagstat_table_collected }

process update_samtools_flagstat_table {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from sambamba_flagstat_table_collected

    output:
    file("${output_file}") into samtools_flagstat_table_ch
    val('') into done_samtools_flagstat_table_update

    script:
    output_file = "flagstat.tsv"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

process sambamba_dedup {
    // deduplicate alignments
    publishDir "${params.outputDir}/alignments/deduplicated", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam) from samples_bam2

    output:
    set val(sampleID), file("${bam_file}") into samples_dd_bam, samples_dd_bam2
    set val(sampleID), file("${bam_file}"), file("${bai_file}") into samples_dd_bam3, samples_dd_bam4, samples_dd_bam5
    set val(sampleID), file("${log_file}") into sambamba_dedup_logs
    val(sampleID) into done_sambamba_dedup

    script:
    prefix = "${sampleID}"
    bam_file = "${prefix}.dd.bam"
    bai_file = "${prefix}.dd.bam.bai"
    log_file = "${prefix}.dd.log"
    """
    sambamba markdup \
    --remove-duplicates \
    --nthreads \${NSLOTS:-\${NTHREADS:-1}} \
    --hash-table-size 525000 \
    --overflow-list-size 525000 \
    "${sample_bam}" "${bam_file}"

    # make a copy of the .command.err Nextflow log file for parsing
    cat .command.err > "${log_file}"

    samtools index "${bam_file}"
    """
}

process sambamba_dedup_log_table {
    // convert the dedup stats the a flat table
    publishDir "${params.outputDir}/alignment-stats", mode: 'copy'

    input:
    set val(sampleID), file(log_file) from sambamba_dedup_logs

    output:
    file("${output_file}") into sambamba_dedup_log_tables
    val(sampleID) into done_sambamba_dedup_log_table

    script:
    prefix = "${sampleID}"
    output_file = "${prefix}.dd.tsv"
    """
    dedup-log2table.R "${log_file}" tmp.tsv
    paste-col.py -i tmp.tsv --header "Sample" -v "${sampleID}" > "${output_file}"
    """
}
sambamba_dedup_log_tables.collectFile(name: ".reads.dedup.tsv", keepHeader: true).set { sambamba_dedup_log_tables_collected }

process update_sambamba_dedup_log_table{
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from sambamba_dedup_log_tables_collected

    output:
    file("${output_file}") into sambamba_dedup_log_table_ch
    val('') into done_update_sambamba_dedup_log_table

    script:
    output_file = "reads.dedup.tsv"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

process samtools_dedup_flagstat {
    // dedup alignment stats
    publishDir "${params.outputDir}/alignment-stats", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam) from samples_dd_bam2

    output:
    set val(sampleID), file("${flagstat}") into dedup_flagstats
    val(sampleID) into done_sambamba_dedup_flagstat

    script:
    prefix = "${sampleID}"
    flagstat = "${prefix}.dd.flagstat.txt"
    """
    samtools flagstat "${sample_bam}" > "${flagstat}"
    """

}

process samtools_dedup_flagstat_table {
    // convert dedup stats to a flat table
    publishDir "${params.outputDir}/alignment-stats", mode: 'copy'

    input:
    set val(sampleID), file(flagstat) from dedup_flagstats

    output:
    file("${output_file}") into sambamba_dedup_flagstat_tables
    val(sampleID) into done_sambamba_dedup_flagstat_table

    script:
    prefix = "${sampleID}"
    output_file = "${prefix}.dd.flagstat.tsv"
    """
    flagstat2table.R "${flagstat}" tmp.tsv
    paste-col.py -i tmp.tsv --header "Sample" -v "${sampleID}" > "${output_file}"
    """
}
sambamba_dedup_flagstat_tables.collectFile(name: ".flagstat.dedup.tsv", keepHeader: true).set { sambamba_dedup_flagstat_tables_collected }

process update_samtools_dedup_flagstat_table {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from sambamba_dedup_flagstat_tables_collected

    output:
    file("${output_file}") into samtools_dedup_flagstat_table_ch
    val('') into done_update_samtools_dedup_flagstat_table

    script:
    output_file = "flagstat.dedup.tsv"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}


// setup downstream Channels
samples_dd_bam.combine(ref_fasta)
            .combine(ref_fai)
            .combine(ref_dict)
            .combine(targets_bed)
            .combine(gatk_1000G_phase1_indels_vcf)
            .combine(gatk_1000G_phase1_indels_vcf_idx)
            .combine(mills_and_1000G_gold_standard_indels_vcf)
            .combine(mills_and_1000G_gold_standard_indels_vcf_idx)
            .combine(dbsnp_ref_vcf)
            .combine(dbsnp_ref_vcf_idx)
            .set { samples_dd_bam_ref_gatk }

// MAIN REALIGNMENT AND RECALIBRATION STEP
process gatk_RealignerTargetCreator {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignment-stats", pattern: "${intervals_file}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(gatk_1000G_phase1_indels_vcf_idx), file(mills_and_1000G_gold_standard_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf_idx), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from samples_dd_bam_ref_gatk

    output:
    set val(sampleID), file("${intervals_file}"), file(sample_bam) into realigned_intervals_tables
    val(sampleID) into done_gatk_RealignerTargetCreator

    script:
    prefix = "${sampleID}"
    intervals_file = "${prefix}.RealignerTargetCreator.intervals"
    """
    gatk.sh -T RealignerTargetCreator \
    -dt NONE \
    --logging_level ERROR \
    -nt \${NSLOTS:-\${NTHREADS:-1}} \
    --reference_sequence "${ref_fasta}" \
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${intervals_file}"
    """
}
realigned_intervals_tables.combine(ref_fasta6)
            .combine(ref_fai6)
            .combine(ref_dict6)
            .combine(targets_bed5)
            .combine(gatk_1000G_phase1_indels_vcf2)
            .combine(gatk_1000G_phase1_indels_vcf_idx2)
            .combine(mills_and_1000G_gold_standard_indels_vcf2)
            .combine(mills_and_1000G_gold_standard_indels_vcf_idx2)
            .combine(dbsnp_ref_vcf6)
            .combine(dbsnp_ref_vcf_idx6)
            .set { realigned_intervals_tables_comb }

process gatk_IndelRealigner {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignments/realigned", pattern: "*.dd.ra.bam*", mode: 'copy'

    input:
    set val(sampleID), file(intervals_file), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(gatk_1000G_phase1_indels_vcf_idx), file(mills_and_1000G_gold_standard_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf_idx), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from realigned_intervals_tables_comb

    output:
    set val(sampleID), file("${ra_bam_file}"), file("${ra_bai_file}") into realigned_intervals_bams
    val(sampleID) into done_gatk_IndelRealigner

    script:
    prefix = "${sampleID}"
    ra_bam_file = "${prefix}.dd.ra.bam"
    ra_bai_file = "${prefix}.dd.ra.bam.bai"
    """
    samtools index "${sample_bam}"

    gatk.sh -T IndelRealigner \
    -dt NONE \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    --maxReadsForRealignment 50000 \
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    -targetIntervals "${intervals_file}" \
    --input_file "${sample_bam}" \
    --out "${ra_bam_file}"

    samtools index "${ra_bam_file}"
    """
}
realigned_intervals_bams.combine(ref_fasta7)
            .combine(ref_fai7)
            .combine(ref_dict7)
            .combine(targets_bed6)
            .combine(gatk_1000G_phase1_indels_vcf3)
            .combine(gatk_1000G_phase1_indels_vcf_idx3)
            .combine(mills_and_1000G_gold_standard_indels_vcf3)
            .combine(mills_and_1000G_gold_standard_indels_vcf_idx3)
            .combine(dbsnp_ref_vcf7)
            .combine(dbsnp_ref_vcf_idx7)
            .set { realigned_intervals_bams_comb }

process gatk_BaseRecalibrator {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignments/stats", pattern: "*.table1.txt", mode: 'copy'

    input:
    set val(sampleID), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(gatk_1000G_phase1_indels_vcf_idx), file(mills_and_1000G_gold_standard_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf_idx), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from realigned_intervals_bams_comb

    output:
    set val(sampleID), file("${table1}"), file(ra_bam_file), file(ra_bai_file) into recalibrated_bases_table1
    val(sampleID) into done_gatk_BaseRecalibrator

    script:
    prefix = "${sampleID}"
    table1 = "${prefix}.BaseRecalibrator.table1.txt"
    """
    gatk.sh -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-\${NTHREADS:-1}} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${ra_bam_file}" \
    --out "${table1}"
    """
}
recalibrated_bases_table1.combine(ref_fasta8)
            .combine(ref_fai8)
            .combine(ref_dict8)
            .tap { recalibrated_bases_table1_ref }
            .combine(targets_bed7)
            .combine(gatk_1000G_phase1_indels_vcf4)
            .combine(gatk_1000G_phase1_indels_vcf_idx4)
            .combine(mills_and_1000G_gold_standard_indels_vcf4)
            .combine(mills_and_1000G_gold_standard_indels_vcf_idx4)
            .combine(dbsnp_ref_vcf8)
            .combine(dbsnp_ref_vcf_idx8)
            .set { recalibrated_bases_table1_comb }

process gatk_BaseRecalibratorBQSR {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignments/stats", pattern: "${table2}", mode: 'copy'

    input:
    set val(sampleID), file(table1), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(gatk_1000G_phase1_indels_vcf_idx), file(mills_and_1000G_gold_standard_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf_idx), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from recalibrated_bases_table1_comb

    output:
    set val(sampleID), file(table1), file("${table2}"), file(ra_bam_file), file(ra_bai_file) into recalibrated_bases_table2
    val(sampleID) into done_gatk_BaseRecalibratorBQSR

    script:
    prefix = "${sampleID}"
    table2 = "${prefix}.BaseRecalibrator.table2.txt"
    """
    gatk.sh -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-\${NTHREADS:-1}} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${ra_bam_file}" \
    -BQSR "${table1}" \
    --out "${table2}"
    """
}
recalibrated_bases_table2.combine(ref_fasta9)
            .combine(ref_fai9)
            .combine(ref_dict9)
            .set { recalibrated_bases_table2_comb }

process gatk_AnalyzeCovariates {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignments/stats", mode: 'copy'
    // sometimes this breaks on small datasets
    // https://gatkforums.broadinstitute.org/gatk/discussion/4213/analyzecovariates-error
    // currently not used in any downstream process so ignore if it fails to run
    errorStrategy = "ignore"

    input:
    set val(sampleID), file(table1), file(table2), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict) from recalibrated_bases_table2_comb

    output:
    file("${csv_file}")
    file("${pdf_file}")

    script:
    prefix = "${sampleID}"
    csv_file = "${prefix}.AnalyzeCovariates.csv"
    pdf_file = "${prefix}.AnalyzeCovariates.pdf"
    """
    gatk.sh -T AnalyzeCovariates \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    -before "${table1}" \
    -after "${table2}" \
    -csv "${csv_file}" \
    -plots "${pdf_file}"
    """
}

process gatk_PrintReads {
    // re-align and recalibrate alignments for later variant calling
    publishDir "${params.outputDir}/alignments/recalibrated", mode: 'copy'

    input:
    set val(sampleID), file(table1), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict) from recalibrated_bases_table1_ref

    output:
    set val(sampleID), file("${ra_rc_bam_file}"), file("${ra_rc_bai_file}") into samples_dd_ra_rc_bam, samples_dd_ra_rc_bam2, samples_dd_ra_rc_bam3, samples_dd_ra_rc_bam4, samples_dd_ra_rc_bam5
    val(sampleID) into done_gatk_PrintReads

    script:
    prefix = "${sampleID}"
    ra_rc_bam_file = "${prefix}.dd.ra.rc.bam"
    ra_rc_bai_file = "${prefix}.dd.ra.rc.bam.bai"
    """
    gatk.sh -T PrintReads \
    --logging_level ERROR \
    -nct \${NSLOTS:-\${NTHREADS:-1}} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -BQSR "${table1}" \
    --input_file "${ra_bam_file}" \
    --out "${ra_rc_bam_file}"

    samtools index "${ra_rc_bam_file}"
    """
}

// setup downstream Channels for per-sample analyses
samples_dd_ra_rc_bam.combine(ref_fasta2)
                    .combine(ref_fai2)
                    .combine(ref_dict2)
                    .tap { samples_dd_ra_rc_bam_ref_nointervals }
                    .combine(targets_bed2)
                    .into { samples_dd_ra_rc_bam_ref;
                            samples_dd_ra_rc_bam_ref2;
                            samples_dd_ra_rc_bam_ref3;
                            samples_dd_ra_rc_bam_ref4;
                            samples_dd_ra_rc_bam_ref5;
                            samples_dd_ra_rc_bam_ref6;
                            samples_dd_ra_rc_bam_ref7;
                            samples_dd_ra_rc_bam_ref8;
                            samples_dd_ra_rc_bam_ref9;
                            samples_dd_ra_rc_bam_ref10 }


samples_dd_ra_rc_bam_ref.combine( dbsnp_ref_vcf3 )
                        .combine(dbsnp_ref_vcf_idx3)
                        .into { samples_dd_ra_rc_bam_ref_dbsnp;
                                samples_dd_ra_rc_bam_ref_dbsnp2 }



process qc_target_reads_gatk_genome {
    // calculate metrics on coverage across whole genome
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_ra_rc_bam_ref_nointervals

    output:
    set val(sampleID), val(mode), file("${sample_summary}") into qc_target_reads_gatk_genomes
    file("${sample_statistics}")
    val(sampleID) into done_qc_target_reads_gatk_genome

    script:
    mode = "genome"
    prefix = "${sampleID}.${mode}"
    sample_statistics = "${prefix}.sample_statistics"
    sample_summary = "${prefix}.sample_summary"
    """
    gatk.sh -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-\${NTHREADS:-1}} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 \
    -ct 50 \
    -ct 100 \
    -ct 200 \
    -ct 300 \
    -ct 400 \
    -ct 500 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${prefix}"
    """
}


process qc_target_reads_gatk_pad500 {
    // calculate metrics on coverage of target regions +/- 500bp
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref2

    output:
    set val(sampleID), val(mode), file("${sample_summary}") into qc_target_reads_gatk_pad500s
    file("${sample_statistics}")
    val(sampleID) into done_qc_target_reads_gatk_pad500

    script:
    mode = "pad500"
    prefix = "${sampleID}.${mode}"
    sample_statistics = "${prefix}.sample_statistics"
    sample_summary = "${prefix}.sample_summary"
    """
    gatk.sh -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 \
    -ct 50 \
    -ct 100 \
    -ct 200 \
    -ct 300 \
    -ct 400 \
    -ct 500 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 500 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${prefix}"
    """
}

process qc_target_reads_gatk_pad100 {
    // calculate metrics on coverage of target regions +/- 100bp
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref3

    output:
    set val(sampleID), val(mode), file("${sample_summary}") into qc_target_reads_gatk_pad100s
    file("${sample_statistics}")
    val(sampleID) into done_qc_target_reads_gatk_pad100

    script:
    mode = "pad100"
    prefix = "${sampleID}.${mode}"
    sample_statistics = "${prefix}.sample_statistics"
    sample_summary = "${prefix}.sample_summary"
    """
    gatk.sh -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 \
    -ct 50 \
    -ct 100 \
    -ct 200 \
    -ct 300 \
    -ct 400 \
    -ct 500 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 100 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${prefix}"
    """
}

process qc_target_reads_gatk_bed {
    // calculate metrics on coverage of target regions; exact region only
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref5

    output:
    set val(sampleID), val(mode), file("${sample_summary}") into qc_target_reads_gatk_beds
    set val(sampleID), val(mode), file("${sample_interval_summary}") into qc_target_reads_gatk_beds_intervals
    val(sampleID) into done_qc_target_reads_gatk_bed

    script:
    mode = "bed"
    prefix = "${sampleID}.${mode}"
    sample_summary = "${prefix}.sample_summary"
    sample_statistics = "${prefix}.sample_statistics"
    sample_interval_summary = "${prefix}.sample_interval_summary"
    sample_interval_statistics = "${prefix}.sample_interval_statistics"
    sample_cumulative_coverage_proportions = "${prefix}.sample_cumulative_coverage_proportions"
    sample_cumulative_coverage_counts = "${prefix}.sample_cumulative_coverage_counts"
    """
    gatk.sh -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    --logging_level ERROR \
    --omitDepthOutputAtEachBase \
    -ct 10 \
    -ct 50 \
    -ct 100 \
    -ct 200 \
    -ct 300 \
    -ct 400 \
    -ct 500 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --nBins 999 \
    --start 1 \
    --stop 1000 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${prefix}"
    """
    // TODO: change this to levels: 10, 50, 100, 200, 300, 400, 500
}

process update_coverage_tables {
    // add more information to the region coverage output
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), val(mode), file(sample_summary) from qc_target_reads_gatk_beds.concat(qc_target_reads_gatk_pad100s, qc_target_reads_gatk_pad500s, qc_target_reads_gatk_genomes)

    output:
    file("${output_file}") into updated_coverage_tables
    val(sampleID) into done_update_coverage_tables

    script:
    prefix = "${sampleID}.${mode}"
    output_file = "${prefix}.sample_summary.tsv"
    """
    sample-summary2table.R "${sample_summary}" tmp.tsv
    paste-col.py -i tmp.tsv --header "Mode" -v "${mode}" > "${output_file}"
    """
}
updated_coverage_tables.collectFile(name: ".${sample_coverage_file}", keepHeader: true).set { updated_coverage_tables_collected }

process update_updated_coverage_tables_collected {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from updated_coverage_tables_collected

    output:
    file("${output_file}") into sample_coverage_file_ch
    val('') into done_update_updated_coverage_tables_collected

    script:
    output_file = "${sample_coverage_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

process coverage_intervals_to_table {
    // convert the intervals coverage file into a formatted table for downstream usage
    input:
    set val(sampleID), val(mode), file(sample_interval_summary) from qc_target_reads_gatk_beds_intervals

    output:
    set val(sampleID), val(mode), file("${output_file}") into updated_coverage_interval_tables

    script:
    prefix = "${sampleID}.${mode}"
    output_file = "${prefix}.sample_interval_summary.table.tsv"
    """
    sample-interval-summary2table.R "${sample_interval_summary}" "${output_file}"
    """
}

cov_ANNOVAR_BUILD_VERSION = "hg19"
cov_ANNOVAR_PROTOCOL = "refGene"
cov_ANNOVAR_OPERATION = "g"
process annotate_coverage_intervals {
    // annotate the coverage table regions
    input:
    set val(sampleID), val(mode), file(interval_table), file(annovar_db_dir) from updated_coverage_interval_tables.combine(annovar_db_dir3)

    output:
    set val(sampleID), val(mode), file("${annotations_tsv}") into annotated_interval_tables

    script:
    prefix = "${sampleID}.${mode}.sample_interval_summary"
    noHeader_tsv = "${prefix}.noHeader.tsv"
    interval_tmp = "${prefix}.intervals.tmp"
    remainder_tsv = "${prefix}.remainder.tsv"
    avinput_file = "${prefix}.avinput"
    // avinput_tsv = "${prefix}.avinput.tsv"
    annovar_output_txt = "${prefix}.${cov_ANNOVAR_BUILD_VERSION}_multianno.txt"
    // annovar_output_vcf = "${prefix}.${cov_ANNOVAR_PROTOCOL}_multianno.vcf"
    annotations_tsv = "${prefix}.annotations.tsv"
    """
    # convert table to ANNOVAR format for annotation; http://annovar.openbioinformatics.org/en/latest/user-guide/input/
    # strip headers and add '0' cols for ref and alt
    tail -n +2 "${interval_table}" > "${noHeader_tsv}"
    cut -f1-3 "${noHeader_tsv}" > "${interval_tmp}"
    cut -f4- "${noHeader_tsv}" > "${remainder_tsv}"
    sed -e 's|\$|\t0|' -i "${interval_tmp}"
    sed -e 's|\$|\t0|' -i "${interval_tmp}"
    paste "${interval_tmp}" "${remainder_tsv}" > "${avinput_file}"

    # annotate
    table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
    --buildver "${cov_ANNOVAR_BUILD_VERSION}" \
    --remove \
    --protocol "${cov_ANNOVAR_PROTOCOL}" \
    --operation "${cov_ANNOVAR_OPERATION}" \
    --nastring . \
    --outfile "${prefix}"

    merge-interval-tables.R "${interval_table}" "${annovar_output_txt}" "${avinput_file}" "${annotations_tsv}"
    """
}

process update_interval_tables {
    // add more information to the intervals coverage output
    publishDir "${params.outputDir}/metrics/coverage/${mode}", mode: 'copy'

    input:
    set val(sampleID), val(mode), file(sample_interval_summary) from annotated_interval_tables

    output:
    file("${output_file}") into updated_annotated_interval_tables
    val(sampleID) into done_update_interval_tables

    script:
    prefix = "${sampleID}.${mode}"
    output_file = "${prefix}.sample_interval_summary.tsv"
    """
    paste-col.py -i "${sample_interval_summary}" --header "Mode" -v "${mode}" | \
    paste-col.py --header "Sample" -v "${sampleID}" > "${output_file}"
    """
}
updated_annotated_interval_tables.collectFile(name: ".${interval_coverage_file}", keepHeader: true).set { updated_coverage_interval_tables_collected }

process update_updated_coverage_interval_tables_collected {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from updated_coverage_interval_tables_collected

    output:
    file("${output_file}") into interval_coverage_file_ch
    val('') into done_update_updated_coverage_interval_tables_collected

    script:
    output_file = "${interval_coverage_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

// not currently used
// process pad_bed {
//     // add 10bp to the regions
//     publishDir "${params.outputDir}/targets", overwrite: true // , mode: 'copy'
//
//     input:
//     set file(targets_bed_file), file(ref_chrom_sizes) from targets_bed3.combine(ref_chrom_sizes)
//
//     output:
//     file("targets.pad10.bed") into targets_pad_bed
//     val("targets.pad10.bed") into done_pad_bed
//
//     script:
//     """
//     cat "${targets_bed_file}" | LC_ALL=C sort -k1,1 -k2,2n | bedtools slop -g "${ref_chrom_sizes}" -b 10 | bedtools merge -d 5 > targets.pad10.bed
//     """
// }

process lofreq {
    // high sensitivity variant calling for low frequency variants
    // NOTE: lofreq_norm_vcfs currently required to start variant channel for rest of pipeline
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref6

    output:
    set val(caller), val(type), val(sampleID), file("${norm_vcf}") into lofreq_norm_vcfs
    file("${vcf_file}")
    file("${norm_vcf}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    val(sampleID) into done_lofreq

    script:
    caller = "LoFreq"
    type = "NA"
    prefix = "${sampleID}.${caller}.${type}"
    vcf_file = "${prefix}.vcf"
    norm_vcf = "${prefix}.norm.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    """
    lofreq call-parallel \
    --call-indels \
    --pp-threads \${NSLOTS:-\${NTHREADS:-1}} \
    --ref "${ref_fasta}" \
    --bed "${targets_bed_file}" \
    --out "${vcf_file}" \
    "${sample_bam}"

    cat ${vcf_file} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """
}

process gatk_hc {
    // variant calling
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"
    // publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    // publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"
    // publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"


    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from samples_dd_ra_rc_bam_ref_dbsnp2

    output:
    set val(caller), val(type), val(sampleID), file("${norm_vcf}") into sample_vcf_hc
    file("${vcf_file}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    file("${norm_vcf}")
    val(sampleID) into done_gatk_hc

    script:
    caller = "HaplotypeCaller"
    type = "NA"
    prefix = "${sampleID}.${caller}.${type}"
    vcf_file = "${prefix}.vcf"
    norm_vcf = "${prefix}.norm.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    """
    gatk.sh -T HaplotypeCaller \
    -dt NONE \
    --logging_level ERROR \
    -nct \${NSLOTS:-\${NTHREADS:-1}} \
    --max_alternate_alleles 3 \
    --standard_min_confidence_threshold_for_calling 50 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${vcf_file}"

    cat ${vcf_file} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """
}

process varscan_snp {
    // https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
    // https://www.biostars.org/p/93673/
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_snp_output}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref9

    output:
    file("${vcf_snp_output}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    set val(caller), val(type), val(sampleID), file("${norm_vcf}") into varscan_snp_vcfs

    when:
    disable_varscan2 != true

    script:
    caller = "VarScan2"
    type = "snp"
    prefix = "${sampleID}.${caller}.${type}"
    vcf_snp_output = "${prefix}.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    norm_vcf = "${prefix}.norm.vcf"
    vcf_sample_list = "vcf-sample-list.txt"
    """
    # make vcf-sample-list
    echo "${sampleID}" > "${vcf_sample_list}"

    # VarScan2 with default settings
    samtools mpileup \
    --no-BAQ \
    --positions "${targets_bed_file}" \
    --fasta-ref "${ref_fasta}" \
    "${sample_bam}" | \
    varscan.sh mpileup2snp \
    --min-coverage 8 \
    --min-reads2 2 \
    --min-avg-qual 15 \
    --min-var-freq 0.01 \
    --min-freq-for-hom 0.75 \
    --p-value 0.99 \
    --strand-filter 1 \
    --vcf-sample-list "${vcf_sample_list}" \
    --output-vcf > "${vcf_snp_output}"

    # normalize and split vcf entries
    cat ${vcf_snp_output} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """

}

process varscan_indel {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_indel_output}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref10

    output:
    file("${vcf_indel_output}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    set val(caller), val(type), val(sampleID), file("${norm_vcf}") into varscan_indel_vcfs

    when:
    disable_varscan2 != true

    script:
    caller = "VarScan2"
    type = "indel"
    prefix = "${sampleID}.${caller}.${type}"
    samtools_mpileup_output = "${prefix}.mpileup"
    vcf_indel_output = "${prefix}.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    norm_vcf = "${prefix}.norm.vcf"
    vcf_sample_list = "vcf-sample-list.txt"
    """
    # make vcf-sample-list
    echo "${sampleID}" > "${vcf_sample_list}"

    # VarScan2 with default settings
    samtools mpileup \
    --no-BAQ \
    --positions "${targets_bed_file}" \
    --fasta-ref "${ref_fasta}" \
    "${sample_bam}" | \
    varscan.sh mpileup2indel \
    --min-coverage 8 \
    --min-reads2 2 \
    --min-avg-qual 15 \
    --min-var-freq 0.01 \
    --min-freq-for-hom 0.75 \
    --p-value 0.99 \
    --strand-filter 1 \
    --vcf-sample-list "${vcf_sample_list}" \
    --output-vcf > "${vcf_indel_output}"

    # normalize and split vcf entries
    cat ${vcf_indel_output} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """
}

// set up channel for filtering unpaired vcfs
lofreq_norm_vcfs.mix(varscan_indel_vcfs, varscan_snp_vcfs, sample_vcf_hc)
    .combine(ref_fasta16)
    .combine(ref_fai16)
    .combine(ref_dict16)
    .set { unpaired_vcfs_ref }
process filter_vcf {
    // filter .vcf files
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"

    input:
    set val(caller), val(type), val(sampleID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from unpaired_vcfs_ref

    output:
    set val(caller), val(type), val(sampleID), file("${filtered_vcf}") into (filtered_vcfs, filtered_vcfs2, filtered_vcfs3) // to tsv

    script:
    prefix = "${sampleID}.${caller}.${type}"
    filtered_vcf = "${prefix}.filtered.vcf"
    if ( caller == "VarScan2" )
        """
        # do not report if frequency is less than 1%
        # automatically filters only PASS variants
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "FREQ > 0.01"  \
        # > "${filtered_vcf}"

        # skip filtering of VarScan2 because it reports the 'FREQ' as a % in the .vcf; "100%", etc.
        # TODO: come up with a filtering method for this

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if ( caller == "LoFreq" )
        """
        # do not report if:
        # - frequency is less than 1%, greater than 99%
        # - depth less than 200
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "AF > 0.01"  \
        # -select "AF < 0.99"  \
        # -select "DP > 100"  \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if ( caller == "HaplotypeCaller" )
        """
        # report if
        # alternate allele freq (allele depth / depth) greater than 0.05 ; 5%
        # more than 5 variant call supporting reads
        # quality reads present (reported depth >0)
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # --sample_name "${sampleID}" \
        # -select "vc.getGenotype('${sampleID}').getAD().1 / vc.getGenotype('${sampleID}').getDP() > 0.05" \
        # -select "vc.getGenotype('${sampleID}').getAD().1 > 5" \
        # -select "vc.getGenotype('${sampleID}').getDP() > 100" \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else
        error "Invalid caller: ${caller}"
}

process vcf_to_tsv {
    // convert .vcf to .tsv
    // TODO: pass .idx file as input here...
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    set val(caller), val(type), val(sampleID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcfs.combine(ref_fasta15).combine(ref_fai15).combine(ref_dict15)

    output:
    set val(caller), val(type), val(sampleID), file(vcf), file("${reformat_tsv}") into (vcf_tsvs, vcf_tsvs2) // to annotation
    set val(caller), val(type), val(sampleID), file("${reformat_tsv}") into vcf_tsvs3


    script:
    prefix = "${sampleID}.${caller}.${type}"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    if ( caller == "VarScan2" )
        """
        # NOTE: The output columns here need to correspond to the columns used in 'reformat-vcf-table.py' !!
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -GF DP \
        -GF AD \
        -GF RD \
        -GF FREQ \
        -GF RBQ \
        -GF ABQ \
        -o "${tsv_file}"

        # reformat and adjust the TSV table for consistency downstream
        reformat-vcf-table.py -c VarScan2 -s "${sampleID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${sampleID}"  | \
        paste-col.py --header "VariantCallerType" -v "${type}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"
        """
        // ##fileformat=VCFv4.1
        // ##source=VarScan2
        // ##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 15">
        // ##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
        // ##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
        // ##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
        // ##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">
        // ##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
        // ##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">
        // ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        // ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        // ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
        // ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
        // ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
        // ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
        // ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
        // ##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
        // ##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
        // ##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
        // ##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
        // ##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
        // ##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
        // ##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1
        // chr1    3102792 .       C       T       .       PASS    ADP=158;WT=0;HET=1;HOM=0;NC=0   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    0/1:0:158:158:156:2:1.27%:9.8E-1:57:56:150:6:2:0
    else if ( caller == "LoFreq" )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F DP \
        -F AF \
        -F SB \
        -F INDEL \
        -F CONSVAR \
        -F HRUN \
        -o "${tsv_file}"

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c LoFreq -s "${sampleID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${sampleID}"  | \
        paste-col.py --header "VariantCallerType" -v "${type}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"
        """
        // ##FILTER=<ID=PASS,Description="All filters passed">
        // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
        // ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
        // ##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
        // ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
        // ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
        // ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
        // ##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
        // ##FILTER=<ID=min_snvqual_72,Description="Minimum SNV Quality (Phred) 72">
        // ##FILTER=<ID=min_indelqual_51,Description="Minimum Indel Quality (Phred) 51">
        // ##FILTER=<ID=min_dp_10,Description="Minimum Coverage 10">
        // ##FILTER=<ID=sb_fdr,Description="Strand-Bias Multiple Testing Correction: fdr corr. pvalue > 0.001000">
        // ##FILTER=<ID=min_snvqual_85,Description="Minimum SNV Quality (Phred) 85">
        // ##FILTER=<ID=min_indelqual_64,Description="Minimum Indel Quality (Phred) 64">
    else if ( caller == "HaplotypeCaller" )
        """
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F FILTER \
        -F QUAL \
        -F AC \
        -F AN \
        -GF AD \
        -GF DP \
        -o "${tsv_file}"

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c HaplotypeCaller -s "${sampleID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${sampleID}" | \
        paste-col.py --header "VariantCallerType" -v "${type}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"
        """
    else
        error "Invalid caller: ${caller}"
}

disable_eval_sample_vcf = true
filtered_vcfs2.combine(dbsnp_ref_vcf4)
    .combine(dbsnp_ref_vcf_idx4)
    .combine(ref_fasta4)
    .combine(ref_fai4)
    .combine(ref_dict4)
    .set { samples_filtered_vcfs }
process eval_sample_vcf {
    // calcaulte variant metrics
    // do we even use this? It takes forever to run and has no downstream dependencies so just disable it for now.
    // only really used for MultQC which is also disabled.
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${eval_file}"

    input:
    set val(caller), val(type), val(sampleID), file(sample_vcf), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(ref_fasta), file(ref_fai), file(ref_dict4) from samples_filtered_vcfs

    output:
    file("${eval_file}")
    val("${sampleID}") into done_eval_sample_vcf

    when:
    disable_eval_sample_vcf != true

    script:
    prefix = "${sampleID}.${caller}.${type}"
    eval_file = "${prefix}.eval.grp"
    """
    gatk.sh -T VariantEval \
    -R "${ref_fasta}" \
    -o "${eval_file}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --eval "${sample_vcf}"
    """
}


// DOWNSTREAM TASKS
//
// DELLY2 SNV STEPS
delly2_modes = [
['deletions', 'DEL'],
['duplications', 'DUP'],
['inversions', 'INV'],
['translocations', 'BND'],
['insertions', 'INS']
]
process delly2 {
    // large structural nucleotide variant calling
    publishDir "${params.outputDir}/snv/${mode}", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref4
    each mode from delly2_modes

    output:
    file "${vcf_file}"
    val(sampleID) into done_delly2

    when:
    disable_delly2 == false

    script:
    type = mode[0]
    arg = mode[1]
    caller = "Delly2"
    prefix = "${sampleID}.${caller}.${type}"
    bcf_file = "${prefix}.bcf"
    vcf_file = "${prefix}.vcf"
    """
    delly call -t "${arg}" -g "${ref_fasta}" -o "${bcf_file}" "${sample_bam}"
    bcftools view "${bcf_file}" > "${vcf_file}"
    """
}


// // REQUIRES ANNOTATIONS FOR DBSNP FILTERING
// TODO: set this back up again
// process vaf_distribution_plot {
//     validExitStatus 0,11 // allow '11' failure triggered by few/no variants
//     errorStrategy 'ignore'
//     publishDir "${params.outputDir}/vaf-distribution-hc", mode: 'copy', overwrite: true
//
//     input:
//     set val(sampleID), file(sample_vcf_annot) from gatk_hc_annotations2
//
//     output:
//     file "${sampleID}_vaf_dist.pdf" into vaf_distribution_plots
//
//     script:
//     """
//     VAF-distribution-plot.R "${sampleID}" "${sample_vcf_annot}"
//     """
// }

// process merge_VAF_plots {
//     executor "local"
//     validExitStatus 0,11 // allow '11' failure triggered by few/no variants
//     errorStrategy 'ignore'
//     publishDir "${params.outputDir}/", mode: 'copy', overwrite: true
//
//     input:
//     file '*' from vaf_distribution_plots.toList()
//
//     output:
//     file "vaf_distributions.pdf"
//
//     script:
//     """
//     if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
//         gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=vaf_distributions.pdf *
//     else
//         exit 11
//     fi
//     """
// }







// SETUP CHANNELS FOR PAIRED TUMOR-NORMAL STEPS

// ~~~~~~~ HAPMAP POOL ANALYSIS ~~~~~~~ //
// get the HapMap sample .bam files into a separate channel
samples_dd_ra_rc_bam5.combine(hapmap_sample_ids).filter { items ->
    def sampleID = items[0]
    def bam = items[1]
    def bai = items[2]
    def HapMapID = items[3]
    sampleID == HapMapID
}.combine(hapmap_pool_ch)
.map { sampleID, bam, bai, HapMapID, HapMapPoolBam, HapMapPoolBai ->
    def HapMapPoolID = "HapMap-Pool"
    def comparisonID = "${sampleID}_${HapMapPoolID}"
    // comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai
    return([ comparisonID, sampleID, bam, bai, HapMapPoolID, HapMapPoolBam, HapMapPoolBai ])
}
.set { hapmap_samples_pool }

// get all the combinations of samples & pairs
samples_dd_ra_rc_bam2.combine(samples_pairs) // [ sampleID, sampleBam, sampleBai, tumorID, normalID ]
                    .filter { item -> // only keep combinations where sample is same as tumor pair sample
                        def sampleID = item[0]
                        def sampleBam = item[1]
                        def sampleBai = item[2]
                        def tumorID = item[3]
                        sampleID == tumorID
                    }
                    .map { item -> // re-order the elements
                        def sampleID = item[0]
                        def sampleBam = item[1]
                        def sampleBai = item[2]
                        def tumorID = item[3]
                        def normalID = item[4]

                        def tumorBam = sampleBam
                        def tumorBai = sampleBai

                        return [ tumorID, tumorBam, tumorBai, normalID ]
                    }
                    .combine(samples_dd_ra_rc_bam3) // combine again to get the samples & files again
                    .filter { item -> // keep only combinations where the normal ID matches the new sample ID
                        def tumorID = item[0]
                        def tumorBam = item[1]
                        def tumorBai = item[2]
                        def normalID = item[3]
                        def sampleID = item[4]
                        def sampleBam = item[5]
                        def sampleBai = item[6]
                        normalID == sampleID
                    }
                    .map {item -> // re arrange the elements
                        def tumorID = item[0]
                        def tumorBam = item[1]
                        def tumorBai = item[2]
                        def normalID = item[3]
                        def sampleID = item[4]
                        def sampleBam = item[5]
                        def sampleBai = item[6]

                        def normalBam = sampleBam
                        def normalBai = sampleBai
                        def comparisonID = "${tumorID}_${normalID}"
                        return [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai ]
                    }
                    .tap { samples_dd_ra_rc_bam_noHapMap_pairs; samples_dd_ra_rc_bam_noHapMap_pairs2 }// added for gatk4
                    .mix(hapmap_samples_pool)
                    .tap { samples_dd_ra_rc_bam_pairs } // make a channel for just this set of data
                    .combine(ref_fasta3) // add reference genome and targets
                    .combine(ref_fai3)
                    .combine(ref_dict3)
                    .tap { samples_dd_ra_rc_bam_pairs_ref_noTargets }
                    .combine(targets_bed4)
                    .tap { samples_dd_ra_rc_bam_pairs_refFasta } // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) ]
                    .tap {  samples_dd_ra_rc_bam_pairs_ref; // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) ]
                            samples_dd_ra_rc_bam_pairs2;
                            samples_dd_ra_rc_bam_pairs_ref3;
                            samples_dd_ra_rc_bam_pairs_ref4;
                            samples_dd_ra_rc_bam_pairs_ref5 // added for gatk4
                          }
                    .combine(microsatellites) // add MSI ref
                    .set { samples_dd_ra_rc_bam_pairs_ref_msi }

samples_dd_bam3.combine(samples_pairs4) // [ sampleID, sampleBam, sampleBai, tumorID, normalID ]
    .filter { items -> // only keep combinations where sample is same as tumor pair sample
        def sampleID = items[0]
        def sampleBam = items[1]
        def sampleBai = items[2]
        def tumorID = items[3]
        sampleID == tumorID
    }
    .map { items -> // re-order the elements
        def sampleID = items[0]
        def sampleBam = items[1]
        def sampleBai = items[2]
        def tumorID = items[3]
        def normalID = items[4]

        def tumorBam = sampleBam
        def tumorBai = sampleBai

        return [ tumorID, tumorBam, tumorBai, normalID ]
    }
    .combine(samples_dd_bam4) // combine again to get the samples & files again
    .filter { items -> // keep only combinations where the normal ID matches the new sample ID
        def tumorID = items[0]
        def tumorBam = items[1]
        def tumorBai = items[2]
        def normalID = items[3]
        def sampleID = items[4]
        def sampleBam = items[5]
        def sampleBai = items[6]
        normalID == sampleID
    }
    .map {items -> // re arrange the elements
        def tumorID = items[0]
        def tumorBam = items[1]
        def tumorBai = items[2]
        def normalID = items[3]
        def sampleID = items[4]
        def sampleBam = items[5]
        def sampleBai = items[6]

        def normalBam = sampleBam
        def normalBai = sampleBai
        def comparisonID = "${tumorID}_${normalID}"
        return [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai ]
    }
    .tap { samples_dd_bam_noHapMap_pairs;
        samples_dd_bam_noHapMap_pairs2 }
    .combine(ref_fasta19) // add reference genome and targets
    .combine(ref_fai19)
    .combine(ref_dict19)
    .tap { samples_dd_bam_noHapMap_pairs_ref;
        samples_dd_bam_noHapMap_pairs_ref2;
        samples_dd_bam_noHapMap_pairs_ref3
      }
        // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, ref_fasta, ref_fai, ref_dict ]

// get the unique chromosomes in the targets bed file
//  for per-chrom paired variant calling
Channel.fromPath( targetsBed )
            .splitCsv(sep: '\t')
            .map{row ->
                row[0]
            }
            .unique()
            .set{ chroms }

process line_chunk {
    // split the .bed file into new files based on desired lines in each file
    input:
    file(bedFile) from targets_bed10

    output:
    file('*') into line_chunk_ch

    script:
    """
    split-bed-lines.py "${bedFile}" "${numTargetSplitLines}"
    """
}
line_chunk_ch.flatten().into {
    line_chunk_ch2;
    line_chunk_ch3 }


// add the reference .vcf
samples_dd_ra_rc_bam_pairs_ref_noTargets.combine(dbsnp_ref_vcf2)
    .combine(dbsnp_ref_vcf_idx2)
    .combine(cosmic_ref_vcf2)
    .combine(cosmic_ref_vcf_idx2)
    .combine( line_chunk_ch2 )
    .map { comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, ref_fasta, ref_fai, ref_dict, dbsnp_vcf, dbsnp_vcf_idx, cosmic_vcf, cosmic_vcf_idx, targets_bed ->
        // get number at the end of the file basename to denote the chunk Label
        def chunkLabel = "${targets_bed.name}".findAll(/\d*$/)[0]

        // re-arrange output order for downstream processes
        return([ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, ref_fasta, ref_fai, ref_dict, targets_bed, dbsnp_vcf, dbsnp_vcf_idx, cosmic_vcf, cosmic_vcf_idx, chunkLabel ])
    }
    .set { samples_dd_ra_rc_bam_pairs_ref_gatk_chrom }


// REQUIRES PAIRED SAMPLES BAM FILES
process msisensor {
    // detection of microsatellite instability regions
    // disable this since it always breaks on small sample datasets
    // this program is buggy;
    // TODO: find a method to pre-filter based on the number of alignments ??
    // validExitStatus 0,139 // allow '139' failure from small dataset; 23039 Segmentation fault      (core dumped)
    errorStrategy 'ignore'
    publishDir "${params.outputDir}/microsatellites", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(microsatellites) from samples_dd_ra_rc_bam_pairs_ref_msi

    output:
    file "${msisensor_output}"
    file "${msisensor_dis}"
    file "${msisensor_germline}"
    file "${msisensor_somatic}"
    val(comparisonID) into done_msisensor

    when:
    disable_msisensor != true

    script:
    prefix = "${comparisonID}"
    msisensor_output = "${prefix}.msisensor"
    msisensor_dis = "${prefix}.msisensor_dis"
    msisensor_germline = "${prefix}.msisensor_germline"
    msisensor_somatic = "${prefix}.msisensor_somatic"
    """
    msisensor msi -d "${microsatellites}" -n "${normalBam}" -t "${tumorBam}" -e "${targets_bed}" -o "${msisensor_output}" -l 1 -q 1 -b \${NSLOTS:-\${NTHREADS:-1}}
    """
}

process mutect2 {
    // paired tumor-normal variant calling
    // NOTE: all '@RG' read groups in .bam file headers should have the same 'SM:' sample ID value!
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"


    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file("targets.bed"), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(cosmic_ref_vcf), file(cosmic_ref_vcf_idx), val(chunkLabel) from samples_dd_ra_rc_bam_pairs_ref_gatk_chrom

    output:
    set val(caller), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${norm_vcf}") into vcfs_mutect2
    set val(caller), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${norm_vcf}") into samples_mutect3
    file("${vcf_file}")
    file("${norm_vcf}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    val(comparisonID) into done_mutect2

    script:
    caller = "MuTect2"
    callerType = "NA"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    vcf_file = "${prefix}.vcf"
    norm_vcf = "${prefix}.norm.vcf"
    filtered_vcf = "${prefix}.filtered.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    """
    # variant calling
    gatk.sh -T MuTect2 \
    -dt NONE \
    --logging_level WARN \
    --standard_min_confidence_threshold_for_calling 30 \
    --max_alt_alleles_in_normal_count 10 \
    --max_alt_allele_in_normal_fraction 0.05 \
    --max_alt_alleles_in_normal_qscore_sum 40 \
    --reference_sequence "${ref_fasta}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --cosmic "${cosmic_ref_vcf}" \
    --intervals "targets.bed" \
    --interval_padding 10 \
    --input_file:tumor "${tumorBam}" \
    --input_file:normal "${normalBam}" \
    --out "${vcf_file}"

    # normalize and split vcf entries
    cat "${vcf_file}" | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """
    // # subset the target regions for the given chromosome
    // # subset_bed.py "${chunkLabel}" "${targets_bed}" > "${bed_subset}"
}

// samples_dd_ra_rc_bam_pairs_ref5.combine(germline_resource_gz) //[ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) ]
//                                .combine(germline_resource_gz_tbi)
//                                .set { samples_dd_ra_rc_bam_pairs_ref_target_germlinevcf }

samples_dd_ra_rc_bam_noHapMap_pairs2.combine(ref_fasta22) // add reference genome and targets //added for gatk4
                                    .combine(ref_fai22)
                                    .combine(ref_dict22)
                                    .combine(targets_bed15)
                                    .combine(germline_resource_gz)
                                    .combine(germline_resource_gz_tbi)
                                    .set { samples_dd_ra_rc_bam_pairs_ref_target_germlinevcf }

process mutect2_gatk4 { //added for gatk4
  // paired tumor-normal variant calling
  // NOTE: all '@RG' read groups in .bam file headers should have the same 'SM:' sample ID value!
  publishDir "${params.outputDir}/variants/MuTect2_GATK4/raw", mode: 'copy', pattern: "*${vcf_file}"

  input:
  set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(germline_resource_gz), file(germline_resource_gz_tbi)from samples_dd_ra_rc_bam_pairs_ref_target_germlinevcf

  output:
  set val(caller), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), file("${vcf_file}") into vcfs_mutect2_gatk4
  file("${vcf_file}")

  script:
  caller = "MuTect2_GATK4"
  prefix = "${comparisonID}.${caller}"
  vcf_file = "${prefix}.vcf"

  """
    gatk --java-options \"-Xms8G -Xmx8G\" Mutect2 \
    --seconds-between-progress-updates 600 \
    --native-pair-hmm-threads 4 \
    --reference "${ref_fasta}" \
    --germline-resource "${germline_resource_gz}" \
    --dont-use-soft-clipped-bases \
    --max-reads-per-alignment-start 100 \
    --intervals "${targets_bed}" \
    --interval-padding 10 \
    --input "${tumorBam}" \
    --input "${normalBam}" \
    --tumor-sample "${tumorID}" \
    --normal-sample "${normalID}" \
    --output "${vcf_file}"
    """
}

process mutect2_vep { //added for Variant Effect Predictor on mutect2 vcf files
  // paired tumor-normal vep calling
  publishDir "${params.outputDir}/variants/MuTect2_GATK4/raw", mode: 'copy', pattern: "*${vep_vcf_file}"

  input:
  set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file(vcf_file),
   file(gnomAD_vcf),file(gnomAD_tbi), file(vep_cache_dir),file(ref_fasta),file(ref_fai) from vcfs_mutect2_gatk4.combine(gnomAD_sites_vcf)
                                                                         .combine(gnomAD_sites_tbi)
                                                                         .combine(vep_cache_dir)
                                                                         .combine(ref_fasta23)
                                                                         .combine(ref_fai23)

  output:
  set val(caller), val(comparisonID), val(tumorID), val(normalID), file("${vep_vcf_file}") into vep_vcfs_mutect2
  file("${vep_vcf_file}")

  script:
  caller = "MuTect2_GATK4"
  prefix = "${comparisonID}.${caller}"
  vep_vcf_file = "${prefix}.vep.vcf"

  """
      vep \
          --fork 4 \
          --species homo_sapiens \
          --offline \
          --everything \
          --shift_hgvs 1 \
          --check_existing \
          --total_length \
          --allele_number \
          --no_escape \
          --refseq \
          --buffer_size 256 \
          --dir "${vep_cache_dir}" \
          --fasta "${ref_fasta}" \
          --input_file "${vcf_file}" \
          --force_overwrite \
          --custom "${gnomAD_vcf},gnomAD,vcf,exact,0,AF_POPMAX,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS" \
          --vcf \
          --output_file "${vep_vcf_file}"
  """
}

process vcf2maf { //convert a VCF into a Mutation Annotation Format (MAF)
  // variants are annotated to all possible gene isoforms
  publishDir "${params.outputDir}/variants/MuTect2_GATK4/raw", mode: 'copy', pattern: "*${maf_file}"

  input:
  set val(caller), val(comparisonID), val(tumorID), val(normalID), file(vep_vcf_file),
   file(ExAC_vcf),file(ExAC_tbi), file(vep_cache_dir),file(ref_fasta),file(ref_fai) from vep_vcfs_mutect2.combine(ExAC_sites_vcf)
                                                                         .combine(ExAC_sites_tbi)
                                                                         .combine(vep_cache_dir)
                                                                         .combine(ref_fasta24)
                                                                         .combine(ref_fai24)

  output:
  set val(caller), val(comparisonID), val(tumorID), val(normalID), file("${maf_file}") into mafs_mutect2
  file("${maf_file}")

  script:
  caller = "MuTect2_GATK4"
  vepPath = "/opt/variant_effect_predictor_96/ensembl-vep-release-96.3/vep"
  retainInfo = "gnomAD_AF_POPMAX,gnomAD_AF_AFR,gnomAD_AF_AMR,gnomAD_AF_ASJ,gnomAD_AF_EAS,gnomAD_AF_FIN,gnomAD_AF_NFE,gnomAD_AF_OTH,gnomAD_AF_SAS"
  prefix = "${comparisonID}.${caller}"
  maf_file = "${prefix}.maf"

  """
      perl /opt/vcf2maf/vcf2maf.pl \
        --species "homo_sapiens" \
        --ncbi-build "GRCh37" \
        --input-vcf "${vcf_file}" \
        --output-maf "${maf_file}" \
        --maf-center "MuTect2" \
        --tumor-id "${tumorID}" \
        --normal-id "${normalID}" \
        --vcf-tumor-id "${tumorID}" \
        --vcf-normal-id "${normalID}" \
        --vep-path "{$vepPath}" \
        --vep-data "${vep_cache_dir}" \
        --ref-fasta "${ref_fasta}" \
        --filter-vcf "${ExAC_vcf}" \
        --buffer-size 265 \
        --max-filter-ac 10 \
        --retain-info "${retainInfo}" \
        --min-hom-vaf 0.7
  """
}

samples_dd_ra_rc_bam_pairs_refFasta.combine(dbsnp_ref_vcf_gz)
.combine(dbsnp_ref_vcf_gz_tbi)
.set { samples_dd_ra_rc_bam_pairs_ref_dbsnp }
// [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(dbsnp_ref_vcf_gz), file(dbsnp_ref_vcf_gz_tbi) ]

process lofreq_somatic {
    // somatic tumor-normal paired variant calling with LoFreq
    publishDir "${params.outputDir}/variants/${caller}/norm", mode: 'copy', pattern: "*${final_snvs_vcf_norm}"
    publishDir "${params.outputDir}/variants/${caller}/norm", mode: 'copy', pattern: "*${final_indels_vcf_norm}"
    publishDir "${params.outputDir}/variants/${caller}/norm", mode: 'copy', pattern: "*${final_minus_dbsnp_snvs_vcf_norm}"
    publishDir "${params.outputDir}/variants/${caller}/norm", mode: 'copy', pattern: "*${final_minus_dbsnp_indels_vcf_norm}"

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(dbsnp_ref_vcf_gz), file(dbsnp_ref_vcf_gz_tbi) from samples_dd_ra_rc_bam_pairs_ref_dbsnp

    output:
    file("${final_snvs_vcf_gz}")
    file("${final_indels_vcf_gz}")
    // file("${final_minus_dbsnp_snvs_vcf_gz}")
    // file("${final_minus_dbsnp_indels_vcf_gz}")
    set val("${caller}"), val("snvs"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${final_snvs_vcf_norm}") into vcfs_lofreq_somatic_snvs_vcf_norm
    set val("${caller}"), val("indels"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${final_indels_vcf_norm}") into vcfs_lofreq_somatic_indels_vcf_norm
    // set val("${caller}"), val("snvs-minus-dbsnp"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${final_minus_dbsnp_snvs_vcf_norm}") into vcfs_lofreq_somatic_snvs_minus_dbsnp_vcf_norm
    // set val("${caller}"), val("indels-minus-dbsnp"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${final_minus_dbsnp_snvs_vcf_norm}") into vcfs_lofreq_somatic_indels_minus_dbsnp_vcf_norm

    script:
    caller = "LoFreqSomatic"
    chunkLabel = "NA"
    prefix = "${comparisonID}.${caller}"
    lofreq_prefix = "${prefix}."

    final_snvs_vcf_gz = "${lofreq_prefix}somatic_final.snvs.vcf.gz"
    final_indels_vcf_gz = "${lofreq_prefix}somatic_final.indels.vcf.gz"
    final_minus_dbsnp_snvs_vcf_gz = "${lofreq_prefix}somatic_final_minus-dbsnp.snvs.vcf.gz"
    final_minus_dbsnp_indels_vcf_gz = "${lofreq_prefix}somatic_final_minus-dbsnp.indels.vcf.gz"

    final_snvs_vcf_norm = "${lofreq_prefix}somatic_final.snvs.norm.vcf"
    final_indels_vcf_norm = "${lofreq_prefix}somatic_final.indels.norm.vcf"
    final_minus_dbsnp_snvs_vcf_norm = "${lofreq_prefix}somatic_final_minus-dbsnp.snvs.norm.vcf"
    final_minus_dbsnp_indels_vcf_norm = "${lofreq_prefix}somatic_final_minus-dbsnp.indels.norm.vcf"
    """
    # paired tumor-normal somatic variant calling with LoFreq
    lofreq somatic \
    --call-indels \
    -n "${normalBam}" \
    -t "${tumorBam}" \
    -f "${ref_fasta}" \
    --threads \${NSLOTS:-\${NTHREADS:-1}} \
    -o "${lofreq_prefix}" \
    -d "${dbsnp_ref_vcf_gz}" \
    -l "${targets_bed}"

    # split multi-allelic entries and left normalize
    zcat ${final_snvs_vcf_gz} | \
    bcftools norm --multiallelics -both --output-type v - | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
    "${final_snvs_vcf_norm}"

    zcat ${final_indels_vcf_gz} | \
    bcftools norm --multiallelics -both --output-type v - | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
    "${final_indels_vcf_norm}"

    zcat ${final_minus_dbsnp_snvs_vcf_gz} | \
    bcftools norm --multiallelics -both --output-type v - | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
    "${final_minus_dbsnp_snvs_vcf_norm}"

    zcat ${final_minus_dbsnp_indels_vcf_gz} | \
    bcftools norm --multiallelics -both --output-type v - | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
    "${final_minus_dbsnp_indels_vcf_norm}"
    """
}

process manta {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bgz), file(targets_tbi) from samples_dd_bam_noHapMap_pairs_ref.combine(targets_zipped)

    output:
    set val("${caller}"), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${candidateSmallIndels_gz}"), file("${candidateSmallIndels_tbi}") into mantaToStrelka
    file("${candidateSV}")
    file("${diploidSV}")
    file("${somaticSV}")
    file("${candidateSmallIndels}")
    // set val("${caller}"), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${candidateSmallIndels}"), file("${candidateSmallIndels_tbi}"), file("${candidateSV}"), file("${candidateSV_tbi}"), file("${diploidSV}"), file("${diploidSV_tbi}"), file("${somaticSV}"), file("${somaticSV_tbi}") into mantaOutput

    script:
    caller = "Manta"
    chunkLabel = "NA"
    callerType = "NA"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    runDir = "${prefix}.Manta"
    candidateSmallIndels = "${prefix}.candidateSmallIndels.vcf"
    candidateSmallIndels_gz = "${prefix}.candidateSmallIndels.vcf.gz"
    candidateSmallIndels_tbi = "${prefix}.candidateSmallIndels.vcf.gz.tbi"
    candidateSV = "${prefix}.candidateSV.vcf"
    candidateSV_gz = "${prefix}.candidateSV.vcf.gz"
    candidateSV_tbi = "${prefix}.candidateSV.vcf.gz.tbi"
    diploidSV = "${prefix}.diploidSV.vcf"
    diploidSV_gz = "${prefix}.diploidSV.vcf.gz"
    diploidSV_tbi = "${prefix}.diploidSV.vcf.gz.tbi"
    somaticSV = "${prefix}.somaticSV.vcf"
    somaticSV_gz = "${prefix}.somaticSV.vcf.gz"
    somaticSV_tbi = "${prefix}.somaticSV.vcf.gz.tbi"
    """
    configManta.py \
    --normalBam "${normalBam}" \
    --tumorBam "${tumorBam}" \
    --referenceFasta "${ref_fasta}" \
    --runDir "${runDir}" \
    --callRegions "${targets_bgz}" \
    --exome

    python ${runDir}/runWorkflow.py \
    -m local \
    -j \${NSLOTS:-\${NTHREADS:-1}}

    # needed for Strelka
    mv ${runDir}/results/variants/candidateSmallIndels.vcf.gz \
    "${candidateSmallIndels_gz}"
    gunzip -c "${candidateSmallIndels_gz}" > "${candidateSmallIndels}"
    mv ${runDir}/results/variants/candidateSmallIndels.vcf.gz.tbi \
    "${candidateSmallIndels_tbi}"

    mv ${runDir}/results/variants/candidateSV.vcf.gz \
    "${candidateSV_gz}"
    gunzip -c "${candidateSV_gz}" > "${candidateSV}"

    mv ${runDir}/results/variants/diploidSV.vcf.gz \
    "${diploidSV_gz}"
    gunzip -c "${diploidSV_gz}" > "${diploidSV}"

    mv ${runDir}/results/variants/somaticSV.vcf.gz \
    "${somaticSV_gz}"
    gunzip -c "${somaticSV_gz}" > "${somaticSV}"
    """
}

samples_dd_bam_noHapMap_pairs.combine(mantaToStrelka)
.filter { items ->
    def comparisonID = items[0]
    def tumorID = items[1]
    def tumorBam = items[2]
    def tumorBai = items[3]
    def normalID = items[4]
    def normalBam = items[5]
    def normalBai = items[6]
    def MantaCaller = items[7]
    def MantaCallerType = items[8]
    def MantaComparisonID = items[9]
    def MantaTumorID = items[10]
    def MantaNormalID = items[11]
    def MantaChunkLabel = items[12]
    def MantaCandidateSmallIndels = items[13]
    def MantaCandidateSmallIndels_tbi = items[14]

    def comparisonID_match = comparisonID == MantaComparisonID
    def tumorID_match = tumorID == MantaTumorID
    def normalID_match = normalID == MantaNormalID
    def all_matches = [ comparisonID_match, tumorID_match, normalID_match ]
    for ( match in all_matches ){
        if ( match == false ){
            return false
        }
    }
    return true
}
.map { items ->
    def comparisonID = items[0]
    def tumorID = items[1]
    def tumorBam = items[2]
    def tumorBai = items[3]
    def normalID = items[4]
    def normalBam = items[5]
    def normalBai = items[6]
    def MantaCaller = items[7]
    def MantaCallerType = items[8]
    def MantaComparisonID = items[9]
    def MantaTumorID = items[10]
    def MantaNormalID = items[11]
    def MantaChunkLabel = items[12]
    def MantaCandidateSmallIndels = items[13]
    def MantaCandidateSmallIndels_tbi = items[14]

    return([comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, MantaCandidateSmallIndels, MantaCandidateSmallIndels_tbi])
}
.combine(ref_fasta20) // add reference genome and targets
.combine(ref_fai20)
.combine(ref_dict20)
.tap { samples_dd_bam_noHapMap_pairs_manta }

process strelka {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(small_indels), file(small_indels_tbi), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bgz), file(targets_tbi) from samples_dd_bam_noHapMap_pairs_manta.combine(targets_zipped2)

    output:
    set val("${caller}"), val("snvs"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${somatic_snvs}") into strelka_snvs
    set val("${caller}"), val("indel"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${somatic_indels}") into strelka_indels

    script:
    caller = "Strelka"
    chunkLabel = "NA"
    callerType = "NA"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    runDir = "${prefix}.Strelka"
    somatic_indels = "${prefix}.somatic.indels.vcf"
    somatic_indels_gz = "${prefix}.somatic.indels.vcf.gz"
    somatic_indels_tbi = "${prefix}.somatic.indels.vcf.gz.tbi"
    somatic_snvs = "${prefix}.somatic.snvs.vcf"
    somatic_snvs_gz = "${prefix}.somatic.snvs.vcf.gz"
    somatic_snvs_tbi = "${prefix}.somatic.snvs.vcf.gz.tbi"
    """
    configureStrelkaSomaticWorkflow.py \
    --normalBam "${normalBam}" \
    --tumorBam "${tumorBam}" \
    --referenceFasta "${ref_fasta}" \
    --indelCandidates "${small_indels}" \
    --runDir ${runDir} \
    --callRegions "${targets_bgz}" \
    --exome

    python ${runDir}/runWorkflow.py \
    -m local \
    -j \${NSLOTS:-\${NTHREADS:-1}}

    mv ${runDir}/results/variants/somatic.indels.vcf.gz \
    "${somatic_indels_gz}"
    gunzip -c "${somatic_indels_gz}" > "${somatic_indels}"

    mv ${runDir}/results/variants/somatic.snvs.vcf.gz \
    "${somatic_snvs_gz}"
    gunzip -c "${somatic_snvs_gz}" > "${somatic_snvs}"
    """
}



// add the split targets bed files
samples_dd_bam_noHapMap_pairs_ref2.combine(line_chunk_ch3)
.map { comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, ref_fasta, ref_fai, ref_dict, targets_bed ->
    // get number at the end of the file basename to denote the chunk Label
    def chunkLabel = "${targets_bed.name}".findAll(/\d*$/)[0]

    return([ chunkLabel, comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, ref_fasta, ref_fai, ref_dict, targets_bed ])
}.set { samples_dd_bam_noHapMap_pairs_targets }

process pindel {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    set val(chunkLabel), val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) from samples_dd_bam_noHapMap_pairs_targets

    when:
    disable_pindel != true

    output:
    set val("${caller}"), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${output_vcf}") into pindel_vcfs
    file("${output_dir}")


    script:
    caller = "Pindel"
    callerType = "indel"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    config_file = "pindel_config.txt"
    output_dir = "${prefix}.pindel_output"
    output_vcf = "${prefix}.vcf"
    insert_size = 500 // 500bp reported by wet lab for sequencing
    """
    # make config file for Pindel
    # labels must be `TUMOR` and `NORMAL` for downstream scripts
    printf "${tumorBam}\t${insert_size}\tTUMOR\n" > "${config_file}"
    printf "${normalBam}\t${insert_size}\tNORMAL\n" >> "${config_file}"

    mkdir "${output_dir}"

    pindel \
    --fasta "${ref_fasta}" \
    --config-file "${config_file}" \
    --output-prefix "${output_dir}/" \
    --number_of_threads \${NSLOTS:-\${NTHREADS:-1}} \
    --max_range_index 2 \
    --include "${targets_bed}"

    pindel2vcf \
    --pindel_output_root "${output_dir}/" \
    --reference "${ref_fasta}" \
    --reference_name hg19 \
    --reference_date 2012_03_15 \
    --gatk_compatible \
    --vcf "${output_vcf}"
    """
    // Variant types reported by Pindel
    // D = deletion
    // SI = short insertion
    // INV = inversion
    // TD = tandem duplication
    // LI = large insertion
    // BP = unassigned breakpoints
}


strelka_snvs.mix(strelka_indels, pindel_vcfs).set{ raw_vcfs_pairs }
process normalize_vcfs_pairs {
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy'

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from raw_vcfs_pairs.combine(ref_fasta21).combine(ref_fai21).combine(ref_dict21)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${norm_vcf}") into norm_vcfs_pairs

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    norm_vcf = "${prefix}.norm.vcf"
    if( caller == 'Pindel' )
        // using fasta-ref with Pindel causes lowercase nucleotides to be embedded in the vcf which breaks annotation downstream
        // Pindel outputs some exact duplicate entries that need to be removed
        """
        # normalize and split vcf entries, remove duplicates
        cat ${vcf} | \
        bcftools norm --multiallelics -both --output-type v - | \
        bcftools norm --rm-dup both --output-type v - | \
        cat -n | sort -k2 -k1n | uniq -f1 | sort -nk1,1 | cut -f2- > \
        "${norm_vcf}"
        """
    else
        """
        cat ${vcf} | \
        bcftools norm --multiallelics -both --output-type v - | \
        bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
        "${norm_vcf}"
        """
}

// get all the paired sample vcfs for downstream processing
vcfs_mutect2.mix(
    vcfs_lofreq_somatic_snvs_vcf_norm,
    vcfs_lofreq_somatic_indels_vcf_norm,
    // vcfs_lofreq_somatic_snvs_minus_dbsnp_vcf_norm,
    // vcfs_lofreq_somatic_indels_minus_dbsnp_vcf_norm,
    norm_vcfs_pairs).set { vcfs_pairs }

process filter_vcf_pairs {
    // filter the .vcf for tumor-normal pairs
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from vcfs_pairs.combine(ref_fasta17).combine(ref_fai17).combine(ref_dict17)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${filtered_vcf}") into filtered_vcf_pairs // to tsv

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    filtered_vcf = "${prefix}.filtered.vcf"
    if( caller == 'MuTect2' )
        """
        # filter VCF
        # report if:
        # only keep 'PASS' entries

        # get the header
        # grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        # grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :

        # old method
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select 'vc.isNotFiltered()' \
        # > "${filtered_vcf}"

        # other criteria:

        # T frequency is more than 3%
        # ( Tumor Allelic depth alt / (Tumor Allelic depth ref + Tumor Allelic depth alt ) )  > 0.03
        # -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) )  > 0.03" \

        # N frequency is less than 5%
        # ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) )  < 0.05
        # -select "(vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) )  < 0.05" \

        # at least 5 variant call supporting reads
        # Tumor Allelic depth alt > 5
        # -select "vc.getGenotype('TUMOR').getAD().1 > 5" \

        # T frequency is sufficiently higher (5x) than N frequency
        # "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
        # ( Tumor Allelic depth alt / ( Tumor Allelic depth ref + Tumor Allelic depth alt ) ) > ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) ) * 5
        # -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) ) > (vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) ) * 5" \

        # vc.getGenotype('TUMOR').getAD().0 ; Tumor Allelic depth ref
        # vc.getGenotype('TUMOR').getAD().1 ; Tumor Allelic depth alt
        # vc.getGenotype('NORMAL').getAD().0 ; Normal Allelic depth ref
        # vc.getGenotype('NORMAL').getAD().1 ; Normal Allelic depth alt

        # example variant format:
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        # CHROM	chr7
        # POS	2946342
        # ID	.
        # REF	A
        # ALT	G
        # QUAL	.
        # FILTER	clustered_events;homologous_mapping_event
        # INFO	ECNT=16;HCNT=3;MAX_ED=65;MIN_ED=5;NLOD=33.41;TLOD=23.04
        # FORMAT	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1
        # TUMOR	0/1:1333,17:0.013:7:10:0.588:0|1:2946342_A_G:40125,535:641:689
        # NORMAL	0/0:137,0:0.00:0:0:.:0|1:2946342_A_G:3959,0:53:80

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if( caller == 'LoFreqSomatic' )
        """
        # do not report if:
        # - frequency is less than 1%, greater than 99%
        # - depth less than 200
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "AF > 0.01"  \
        # -select "AF < 0.99"  \
        # -select "DP > 100"  \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if( caller == 'Strelka' )
        """
        # only keep 'PASS' entries
        # filter out TQSS_NT=2 https://github.com/Illumina/strelka/issues/65
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "TQSS_NT != 2"  \
        # --excludeFiltered \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
        // ##FILTER=<ID=PASS,Description="All filters passed">
        // ##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
        // ##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
        // ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
        // ##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
        // ##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
        // ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
        // ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
        // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples">
        // ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
        // ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
        // ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor">
        // ##INFO=<ID=SNVSB,Number=1,Type=Float,Description="Somatic SNV site strand bias">
        // ##INFO=<ID=PNOISE,Number=1,Type=Float,Description="Fraction of panel containing non-reference noise at this site">
        // ##INFO=<ID=PNOISE2,Number=1,Type=Float,Description="Fraction of panel containing more than one non-reference noise obs at this site">
        // ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
        // ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
        // ##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
        // ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
        // ##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
        // ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
        // ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
        // ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
        // ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
        // ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
        // ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">
    else if( caller == 'Pindel' )
        """
        # only keep 'PASS' entries
        # only keep SNV less than 500bp long
        # variant must have:
        # Tumor Alt >= 5 reads
        # Tumor Ref >= 100 reads
        # Normal Alt < 2 reads
        # Normal Ref >= 100 reads

        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select 'SVLEN < 500' \
        # -select 'vc.getGenotype("TUMOR").getAD().1 > 4' \
        # -select 'vc.getGenotype("TUMOR").getAD().0 > 99' \
        # -select 'vc.getGenotype("NORMAL").getAD().1 < 2' \
        # -select 'vc.getGenotype("NORMAL").getAD().0 > 99' \
        # > "${filtered_vcf}"

        cp "${vcf}" "${filtered_vcf}"
        """
        // ##FILTER=<ID=PASS,Description="All filters passed">
        // ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
        // ##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
        // ##INFO=<ID=PF,Number=1,Type=Integer,Description="The number of samples carry the variant">
        // ##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
        // ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
        // ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
        // ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
        // ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
        // ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        // ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference depth, how many reads support the reference">
        // ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">
    else
        error "Invalid caller: ${caller}"
}

process vcf_to_tsv_pairs {
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcf_pairs.combine(ref_fasta18).combine(ref_fai18).combine(ref_dict18)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file("${reformat_tsv}") into vcf_tsv_pairs // to annotation
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${reformat_tsv}") into vcf_tsv_pairs2

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    if( caller == 'MuTect2' )
        """
        # convert VCF to TSV
        # NOTE: automatically filters for only PASS entries
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
        -GF AD -GF DP -GF AF \
        -o "${tsv_file}"

        # .vcf field descriptions:
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##INFO=<ID=NLOD,Number=1,Type=String,Description="Normal LOD score">
        ##INFO=<ID=TLOD,Number=1,Type=String,Description="Tumor LOD score">

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}"  | \
        paste-col.py --header "Tumor" -v "${tumorID}"  | \
        paste-col.py --header "Normal" -v "${normalID}"  | \
        paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"

        # make sure that the input and output tables have the same number of rows
        if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
        """
    else if( caller == 'LoFreqSomatic' )
        """
        # convert to tsv format
        # NOTE: automatically filters for only PASS entries
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F DP \
        -F DP4 \
        -F AF \
        -F SB \
        -F UQ \
        -F CONSVAR \
        -F AN \
        -F AC \
        -F HRUN \
        -F INDEL \
        -F UNIQ \
        -F SOMATIC \
        -o "${tsv_file}"

        # .vcf field descriptions:
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
        ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
        ##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
        ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
        ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
        ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">
        ##INFO=<ID=UNIQ,Number=0,Type=Flag,Description="Unique, i.e. not detectable in paired sample">
        ##INFO=<ID=UQ,Number=1,Type=Integer,Description="Phred-scaled uniq score at this position">

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c LoFreq -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}" | \
        paste-col.py --header "Tumor" -v "${tumorID}" | \
        paste-col.py --header "Normal" -v "${normalID}" | \
        paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"

        # make sure that the input and output tables have the same number of rows
        if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
        """
    else if( caller == 'Strelka' )
        if ( callerType == "snvs" )
            """
            # convert to tsv format
            # NOTE: automatically filters for only PASS entries
            gatk.sh -T VariantsToTable \
            -R "${ref_fasta}" \
            -V "${vcf}" \
            -F CHROM \
            -F POS \
            -F ID \
            -F REF \
            -F ALT \
            -F FILTER \
            -F DP \
            -F SOMATIC \
            -F QSS \
            -F MQ \
            -F SNVSB \
            -F SomaticEVS \
            -GF DP \
            -GF AU \
            -GF TU \
            -GF CU \
            -GF GU \
            -o "${tsv_file}"

            # vcf header example for Strelka SNVs
            ##fileformat=VCFv4.1
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##source=strelka
            ##source_version=2.9.10
            ##content=strelka somatic snv calls
            ##priorSomaticSnvRate=0.0001
            ##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
            ##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
            ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
            ##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
            ##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
            ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
            ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples">
            ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
            ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
            ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor">
            ##INFO=<ID=SNVSB,Number=1,Type=Float,Description="Somatic SNV site strand bias">
            ##INFO=<ID=PNOISE,Number=1,Type=Float,Description="Fraction of panel containing non-reference noise at this site">
            ##INFO=<ID=PNOISE2,Number=1,Type=Float,Description="Fraction of panel containing more than one non-reference noise obs at this site">
            ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
            ##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
            ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
            ##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
            ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
            ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
            ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
            ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
            ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
            ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">

            # reformat and adjust the TSV table for consistency downstream
            # add extra columns to the VCF TSV file for downstream
            reformat-vcf-table.py -c StrelkaSomaticSNV -s "${tumorID}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumorID}" | \
            paste-col.py --header "Tumor" -v "${tumorID}" | \
            paste-col.py --header "Normal" -v "${normalID}" | \
            paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
            paste-col.py --header "VariantCaller" -v "${caller}" > \
            "${reformat_tsv}"

            # make sure that the input and output tables have the same number of rows
            tsv_lines="\$(wc -l < "${tsv_file}" )"
            reformat_lines="\$( wc -l < "${reformat_tsv}" )"
            if [ "\$tsv_lines" -ne "\$reformat_lines" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi

            vcf_lines="\$(grep -v '^##' ${vcf} | wc -l)"
            if [ "\$vcf_lines" -ne "\$reformat_lines" ]; then echo "ERROR: reformat table has different number of entries than vcf file!"; exit 1; fi
            """
        else if( callerType == 'indel' )
            """
            # convert to tsv format
            # NOTE: automatically filters for only PASS entries
            gatk.sh -T VariantsToTable \
            -R "${ref_fasta}" \
            -V "${vcf}" \
            -F CHROM \
            -F POS \
            -F ID \
            -F REF \
            -F ALT \
            -F FILTER \
            -F SOMATIC \
            -F MQ \
            -F SomaticEVS \
            -F QSI \
            -GF DP \
            -GF TAR \
            -GF TIR \
            -GF TOR \
            -o "${tsv_file}"

            # vcf header example for Strelka indels
            ##fileformat=VCFv4.1
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##source=strelka
            ##source_version=2.9.10
            ##content=strelka somatic indel calls
            ##priorSomaticIndelRate=1e-06
            ##INFO=<ID=QSI,Number=1,Type=Integer,Description="Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal">
            ##INFO=<ID=TQSI,Number=1,Type=Integer,Description="Data tier used to compute QSI">
            ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
            ##INFO=<ID=QSI_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
            ##INFO=<ID=TQSI_NT,Number=1,Type=Integer,Description="Data tier used to compute QSI_NT">
            ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
            ##INFO=<ID=RU,Number=1,Type=String,Description="Smallest repeating sequence unit in inserted or deleted sequence">
            ##INFO=<ID=RC,Number=1,Type=Integer,Description="Number of times RU repeats in the reference allele">
            ##INFO=<ID=IC,Number=1,Type=Integer,Description="Number of times RU repeats in the indel allele">
            ##INFO=<ID=IHP,Number=1,Type=Integer,Description="Largest reference interrupted homopolymer length intersecting with the indel">
            ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
            ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
            ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
            ##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description="Somatic indel possibly overlaps a second indel.">
            ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1">
            ##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Read depth for tier2">
            ##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
            ##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
            ##FORMAT=<ID=TOR,Number=2,Type=Integer,Description="Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2">
            ##FORMAT=<ID=DP50,Number=1,Type=Float,Description="Average tier1 read depth within 50 bases">
            ##FORMAT=<ID=FDP50,Number=1,Type=Float,Description="Average tier1 number of basecalls filtered from original read depth within 50 bases">
            ##FORMAT=<ID=SUBDP50,Number=1,Type=Float,Description="Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases">
            ##FORMAT=<ID=BCN50,Number=1,Type=Float,Description="Fraction of filtered reads within 50 bases of the indel.">
            ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
            ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">

            reformat-vcf-table.py -c StrelkaSomaticIndel -s "${tumorID}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumorID}" | \
            paste-col.py --header "Tumor" -v "${tumorID}" | \
            paste-col.py --header "Normal" -v "${normalID}" | \
            paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
            paste-col.py --header "VariantCaller" -v "${caller}" > \
            "${reformat_tsv}"

            # make sure that the input and output tables have the same number of rows
            if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
            """
        else
            error "Invalid Strelka callerType: ${callerType}"
    else if( caller == 'Pindel' )
        """
        # convert VCF to TSV
        # NOTE: automatically filters for only PASS entries
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F FILTER \
        -F QUAL \
        -F END \
        -F HOMLEN \
        -F SVLEN \
        -F SVTYPE \
        -GF AD \
        -o "${tsv_file}"

        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
        ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference depth, how many reads support the reference">
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
        ##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
        ##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
        ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
        ##INFO=<ID=PF,Number=1,Type=Integer,Description="The number of samples carry the variant">
        ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
        ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c Pindel -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}"  | \
        paste-col.py --header "Tumor" -v "${tumorID}"  | \
        paste-col.py --header "Normal" -v "${normalID}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"

        # make sure that the input and output tables have the same number of rows
        if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
        """
    else
        error "Invalid caller: ${caller}"
}


samples_mutect3.combine(dbsnp_ref_vcf5)
                .combine(dbsnp_ref_vcf_idx5)
                .combine(ref_fasta5)
                .combine(ref_fai5)
                .combine(ref_dict5)
                .set { pairs_filtered_vcfs }

process eval_pair_vcf {
    // variant quality metrics
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${eval_file}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(pairs_vcf), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(ref_fasta), file(ref_fai), file(ref_dict4) from pairs_filtered_vcfs

    output:
    file("${eval_file}")
    val("${comparisonID}") into done_eval_pair_vcf

    when:
    disable_eval_pair_vcf != true

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    eval_file = "${prefix}.eval.grp"
    """
    gatk.sh -T VariantEval \
    -R "${ref_fasta}" \
    -o "${eval_file}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --eval "${pairs_vcf}"
    """
}

// ~~~~~~ VARIANT ANNOTATION ~~~~~~ //
samples_vcfs_tsvs_good = Channel.create()
samples_vcfs_tsvs_bad = Channel.create()
pairs_vcfs_tsvs_good = Channel.create()
pairs_vcfs_tsvs_bad = Channel.create()

// filter out samples with empty variant tables
vcf_tsvs.choice( samples_vcfs_tsvs_good, samples_vcfs_tsvs_bad ){ items ->
    def caller =  items[0]
    def type = items[1]
    def sampleID =  items[2]
    def sample_vcf =  items[3]
    def sample_tsv =  items[4]
    def output = 1 // bad by default
    long count = Files.lines(sample_tsv).count()
    if (count > 1) output = 0 // good if has >1 lines
    return(output)
}

samples_vcfs_tsvs_bad.map { caller, type, sampleID, sample_vcf, sample_tsv ->
    def reason = "Too few lines in sample_tsv, skipping annotation"
    def output = [sampleID, caller, type, reason, "${sample_vcf},${sample_tsv}"].join('\t')
    return(output)
}.set { samples_vcfs_tsvs_bad_logs }

vcf_tsv_pairs.choice( pairs_vcfs_tsvs_good, pairs_vcfs_tsvs_bad ){ items ->
    def caller = items[0]
    def callerType = items[1]
    def comparisonID = items[2]
    def tumorID = items[3]
    def normalID = items[4]
    def chunkLabel = items[5]
    def sample_vcf = items[6]
    def sample_tsv = items[7]
    def output = 1 // bad by default
    long count = Files.lines(sample_tsv).count()
    if (count > 1) output = 0 // good if has >1 lines
    return(output)
}

pairs_vcfs_tsvs_bad.map { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, sample_vcf, sample_tsv ->
    def reason = "Too few lines in sample_tsv, skipping annotation"
    def output = [comparisonID, tumorID, normalID, chunkLabel, caller, callerType, reason, "${sample_vcf},${sample_tsv}"].join('\t')
    return(output)
}.set { pairs_vcfs_tsvs_bad_logs }

process annotate {
    // annotate variants
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/annotations/${caller}", mode: 'copy', pattern: "*${annotations_tsv}"

    input:
    set val(caller), val(type), val(sampleID), file(sample_vcf), file(sample_tsv), file(annovar_db_dir) from samples_vcfs_tsvs_good.combine(annovar_db_dir)

    output:
    file("${annotations_tsv}") into annotations_tables
    set val(sampleID), val(caller), val(type), file("${annotations_tsv}") into annotations_annovar_tables
    set val(caller), val(type), val(sampleID), file("${annotations_tsv}") into (annotations_annovar_tables2, annotations_annovar_tables3, annotations_annovar_tables4, annotations_annovar_tables5)
    file("${avinput_file}")
    file("${avinput_tsv}")
    val(sampleID) into done_annotate

    script:
    prefix = "${sampleID}.${caller}.${type}"
    avinput_file = "${prefix}.avinput"
    avinput_tsv = "${prefix}.avinput.tsv"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annovar_output_vcf = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    annotations_tsv = "${prefix}.annotations.tsv"
    if( caller == 'HaplotypeCaller' )
        """
        # convert to ANNOVAR format
        convert2annovar.pl \
        -includeinfo \
        -format vcf4 \
        "${sample_vcf}" > \
        "${avinput_file}"

        # annovate
        table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --outfile "${prefix}"

        # add headers to the avinput, just the first columns
        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-10 "${avinput_file}" >>  "${avinput_tsv}"

        # merge the tables together
        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else if( caller == 'LoFreq' )
        """
        # convert to ANNOVAR format
        convert2annovar.pl \
        -includeinfo \
        -format vcf4 \
        "${sample_vcf}" > \
        "${avinput_file}"

        table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --outfile "${prefix}"

        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-10 ${avinput_file} >>  "${avinput_tsv}"

        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else if( caller == 'VarScan2' )
        """
        # convert to ANNOVAR format
        convert2annovar.pl \
        -includeinfo \
        -format vcf4 \
        "${sample_vcf}" > \
        "${avinput_file}"

        table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --outfile "${prefix}"

        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-10 ${avinput_file} >>  "${avinput_tsv}"

        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else
        error "Invalid caller: ${caller}"
}

process annotate_pairs {
    // annotate variants
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/annotations/${caller}", mode: 'copy', pattern: "*${annotations_tsv}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(sample_vcf), file(sample_tsv), file(annovar_db_dir) from pairs_vcfs_tsvs_good.combine(annovar_db_dir2)

    output:
    file("${annotations_tsv}") into annotations_tables_pairs
    val(comparisonID) into done_annotate_pairs

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    avinput_file = "${prefix}.avinput"
    avinput_tsv = "${prefix}.avinput.tsv"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annotations_tsv = "${prefix}.annotations.tsv"
    vcf_gt_mod = "${sample_vcf}.GTmod.vcf" // only used for Strelka
    if( caller == 'MuTect2' )
        """
        # annotate .vcf
        table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --outfile "${prefix}"

        # get values from .avinput file
        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"

        # merge the tables together
        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else if( caller == "LoFreqSomatic" )
        """
        # annotate .vcf
        table_annovar.pl "${sample_vcf}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --outfile "${prefix}"

        # get values from .avinput file
        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"

        # merge the tables together
        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else if( caller == "Strelka" )
        """
        # need to add the GT field to Strelka .vcf files
        vcf-GT-mod.py "${sample_vcf}" "${vcf_gt_mod}"

        # annotate .vcf
        table_annovar.pl "${vcf_gt_mod}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --vcfinput \
        --outfile "${prefix}"

        # get values from .avinput file
        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"

        # merge the tables together
        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
        """
    else if( caller == "Pindel" )
        """
        # Pindel produces a lot of near-duplicate entries, need to take extra steps to merge tables successfully
        # merge errors should be detected by mismatched number of table rows after merge, in R script

        # convert to ANNOVAR format
        convert2annovar.pl \
        -includeinfo \
        -format vcf4old \
        "${sample_vcf}" > \
        "${avinput_file}"

        # Annotate
        table_annovar.pl "${avinput_file}" "${annovar_db_dir}" \
        --buildver "${params.ANNOVAR_BUILD_VERSION}" \
        --remove \
        --protocol "${params.ANNOVAR_PROTOCOL}" \
        --operation "${params.ANNOVAR_OPERATION}" \
        --nastring . \
        --outfile "${prefix}"

        # need to re-associate the original vcf columns with the vcf tsv for merge
        grep -v '^##' "${sample_vcf}" | cut -f8-11 | paste "${sample_tsv}" /dev/stdin > tmp.tsv
        grep -v '^##' "${sample_vcf}" | cut -f8-11 | paste ${annovar_output_txt} /dev/stdin > tmp.hg19_multianno.txt

        # merge the tables together
        merge-vcf-tables-Pindel.R \
        "${sample_vcf}" \
        "tmp.tsv" \
        "tmp.hg19_multianno.txt" \
        "${avinput_file}" \
        "${annotations_tsv}"
        """
    else
        error "Invalid caller: ${caller}"
}

process collect_annotation_tables {
    // combine all variants into a single table

    input:
    file('t*') from annotations_tables.concat(annotations_tables_pairs).collect()

    output:
    file("${output_file}") into collected_annotation_tables
    val('') into done_collect_annotation_tables

    script:
    output_file = "all_annotations.tmp.tsv"
    """
    concat-tables.py t* > "${output_file}"
    """
}

process filter_annotation_table {
    // apply extra filter criteria here which couldnt be easily applied earlier in the pipeline
    input:
    file(tsv) from collected_annotation_tables

    output:
    file("${output_file}") into filtered_annotation_table

    script:
    output_file = "annotations.filtered.tsv"
    """
    filter-annotation-table.py -i "${tsv}" -o "${output_file}"
    """
}

process update_collect_annotation_tables {
    // add workflow labels to the annotation table
    // publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from filtered_annotation_table

    output:
    file("${output_file}") into (all_annotations_file_ch, all_annotations_file_ch2, all_annotations_file_ch3, all_annotations_file_ch4)
    val('') into done_update_collect_annotation_tables

    script:
    output_file = "${all_annotations_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

process split_annotation_table_caller {
    // create a separate annotation table file for each variant caller
    publishDir "${params.outputDir}/annotations", mode: 'copy'

    input:
    file(".annotations.tsv") from all_annotations_file_ch2

    output:
    file("*")

    script:
    """
    split-annotation-table-callers.py ".annotations.tsv"
    """
}

process extract_hapmap_pool_annotations {
    // make a separate table just for HapMap-Pool variants; for output convenience
    publishDir "${params.outputDir}/annotations", mode: 'copy'

    input:
    file(tsv) from all_annotations_file_ch3

    output:
    file("${output_file}") into hapmap_pool_annotations

    script:
    output_file = "${all_HapMapPool_annotations_file}"
    """
    head -1 "${tsv}" > "${output_file}"
    grep -i 'HapMap-Pool' "${tsv}" >> "${output_file}" || :
    """
}

process split_annotation_table_paired {
    // Split the annotations table into separate files
    // based on variant callers
    // and whether the variant caller is designated at being used for
    // "paired" or "unpaired" variant calling (tumor normal pairs)
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file("annotations.tsv") from all_annotations_file_ch4

    output:
    file("${output_paired}")
    file("${output_unpaired}")

    script:
    output_paired = "annotations.paired.tsv"
    output_unpaired = "annotations.unpaired.tsv"
    """
    split-annotation-table-pairs.py "annotations.tsv" "${output_paired}.tmp" "${output_unpaired}.tmp"

    annotation-table-drop-cols.py -i "${output_paired}.tmp" -o "${output_paired}" --type paired
    annotation-table-drop-cols.py -i "${output_unpaired}.tmp" -o "${output_unpaired}" --type unpaired

    # make sure that the output table number of records matches the input
    num_original_annotations="\$(( \$(wc -l < annotations.tsv) -1 ))"
    num_paired_annotations="\$(( \$(wc -l < "${output_paired}") -1 ))"
    num_unpaired_annotations="\$(( \$(wc -l < "${output_unpaired}") -1 ))"
    num_total_annotations="\$(( \$num_paired_annotations + \$num_unpaired_annotations ))"
    if [ \$num_total_annotations -ne \$num_original_annotations ]; then echo "ERROR: annotation table has different number of rows!"; exit 1; fi
    """
}

// ~~~~ Tumor Burden Analysis ~~~~ //
samples_dd_ra_rc_bam4.combine(ref_fasta10)
.combine(ref_fai10)
.combine(ref_dict10)
.combine(targets_bed8)
.set { samples_bam_ref }

process gatk_CallableLoci {
    publishDir "${params.outputDir}/loci", mode: 'copy'

    input:
    set val(sampleID), file(bamfile), file(baifile), file(genomeFa), file(genomeFai), file(genomeDict), file(targets_bed) from samples_bam_ref

    output:
    set val("${sampleID}"), file("${output_summary}") into called_loci
    file("${output_summary}")
    file("${output_bed}")
    file("${output_bed_pass}")

    script:
    prefix = "${sampleID}"
    output_summary = "${prefix}.CallableLoci.summary.txt"
    output_bed = "${prefix}.CallableLoci.bed"
    output_bed_pass = "${prefix}.CallableLoci.pass.bed"
    """
    # minDepth 500 = NGS580 call threshold
    # minMappingQuality, minBaseQuality 20 = NGS580 DepthOfCoverage threshold
    gatk.sh \
    -T CallableLoci \
    -R "${genomeFa}" \
    -I "${bamfile}" \
    -summary "${output_summary}" \
    --minMappingQuality 20 \
    --minBaseQuality 20 \
    --minDepth 500 \
    --intervals "${targets_bed}" \
    -o "${output_bed}"

    grep -E 'CALLABLE|PASS' "${output_bed}" > "${output_bed_pass}" || touch "${output_bed_pass}" # exit code 1 if no matches
    """
}

process callable_loci_table {
    publishDir "${params.outputDir}/loci", mode: 'copy'

    input:
    set val(sampleID), file(summary) from called_loci

    output:
    file("${output_summary}") into loci_tables
    set val(sampleID), file("${output_summary}"), file("${output_txt}") into loci_tables2

    script:
    prefix = "${sampleID}"
    output_summary = "${prefix}.CallableLoci.summary.tsv"
    output_txt = "${prefix}.CallableLoci.txt"
    tmpFile = "tmp"
    """
    callable-loci-table.py "${summary}" "${tmpFile}"
    paste-col.py -i "${tmpFile}" --header "Sample" -v "${sampleID}" > "${output_summary}"
    grep 'CALLABLE' "${tmpFile}" | cut -f2 > "${output_txt}"
    """
}

// only keep files with at least 1 variant for TMB analysis
annotations_annovar_tables.filter { sampleID, caller, type, anno_tsv ->
    def count = anno_tsv.readLines().size()
    count > 1
}.filter { sampleID, caller, type, anno_tsv ->
    // dont use indels
    type != "indel"
}.set { annotations_annovar_tables_filtered }

process tmb_filter_variants {
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/tmb/${caller}/annotations", mode: 'copy', pattern: "*${output_variants}"

    input:
    set val(sampleID), val(caller), val(type), file(anno_tsv) from annotations_annovar_tables_filtered

    output:
    set val(sampleID), val(caller), val(type), file("${output_variants}") into tmb_filtered_variants2

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_variants = "${prefix}.annotations.tmb.filtered.tsv"
    """
    tmb-variant-filter.py -c "${caller}" -i "${anno_tsv}" | \
    paste-col.py --header "Sample" -v "${sampleID}" | \
    paste-col.py --header "VariantCallerType" -v "${type}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${output_variants}"
    """
}

loci_tables2.map{ sampleID, summary_table, callable_txt ->
    // read the number of callable loci from the file, just pass the value
    def callable_loci = callable_txt.readLines()[0]
    return([sampleID, summary_table, "${callable_loci}"])
}
// combine against the filtered tmb variants
.combine(tmb_filtered_variants2)
.filter{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, type, filtered_annotations ->
    // only keep the matches between channels
    sampleID_anno == sampleID_loci
}
.map{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, type, filtered_annotations ->
    // read the number of variants from the file
    def num_variants = filtered_annotations.readLines().size() - 1
    return([sampleID_loci, caller, type, callable_loci, num_variants])
}
.filter{ sampleID_loci, caller, type, callable_loci, num_variants ->
    // only keep samples that had callable loci
    callable_loci > 0
}
.set { loci_annotations }

process calculate_tmb {
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/tmb/${caller}/values", mode: 'copy', pattern: "*${output_tmb}"

    input:
    set val(sampleID), val(caller), val(type), val(loci), val(variants) from loci_annotations

    output:
    file("${output_tmb}") into tmbs

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_tmb = "${prefix}.tmb.tsv"
    """
    tmb=\$( calc-tmb.py ${variants} ${loci} )
    printf 'SampleID\tVariantCaller\tVariantCallerType\tnBases\tnVariants\tTMB\n' > "${output_tmb}"
    printf "${sampleID}\t${caller}\t${type}\t${loci}\t${variants}\t\${tmb}\n" >> "${output_tmb}"
    """
}
tmbs.collectFile(name: "${tmb_file}", keepHeader: true, storeDir: "${params.outputDir}").set { tmbs_collected }








// Genomic Signatures
process signatures_variant_filter {
    // need to apply some filter criteria to the variant tables before using them to generate signatures
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/signatures/${caller}/annotations", mode: 'copy', pattern: "*${output_file}"

    input:
    set val(caller), val(type), val(sampleID), file(tsv) from annotations_annovar_tables2

    output:
    set val(caller), val(type), val(sampleID), file("${output_file}") into filtered_signatures_tsvs

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_file = "${prefix}.filtered.tsv"
    if ( caller == "VarScan2" )
        """
        signatures-variant-filter.py -c VarScan2 -i "${tsv}" -o "${output_file}"
        """
    else if ( caller == "HaplotypeCaller" )
        """
        signatures-variant-filter.py -c HaplotypeCaller -i "${tsv}" -o "${output_file}"
        """
    else if ( caller == "LoFreq" )
        """
        signatures-variant-filter.py -c LoFreq -i "${tsv}" -o "${output_file}"
        """
    else
        error "Invalid caller: ${caller}"
}

// set up channels for filtering
sample_sig_good = Channel.create()
sample_sig_bad = Channel.create()

// minimum number of entries required to run signatures
deconstructSigs_variant_min = 55

filtered_signatures_tsvs.choice( sample_sig_good, sample_sig_bad ){ items ->
    // make sure there are enough variants in the VCF to proceed!
    def caller = items[0]
    def type = items[1]
    def sampleID = items[2]
    def tsv = items[3]
    def output_ch = 1 // sample_sig_bad
    long count = Files.lines(tsv).count()
    if ( count > deconstructSigs_variant_min ) {
        output_ch = 0 // sample_sig_good
    }
    return(output_ch)
}

sample_sig_bad.map {  caller, type, sampleID, filtered_vcf ->
    def reason = "Fewer than ${deconstructSigs_variant_min} variants in .tsv file, skipping genomic signatures"
    def output = [sampleID, caller, type, reason, filtered_vcf].join('\t')
    return(output)
}.set { sample_sig_bad_logs }

// dont make signatures for indels
sample_sig_good.filter { items ->
    def caller = items[0]
    def type = items[1]
    def sampleID = items[2]
    def tsv = items[3]

    type != 'indel'
}.set { sample_sig_good_filtered }

process deconstructSigs_signatures {
    // search for mutation signatures
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/signatures/${caller}/Rds", mode: 'copy', pattern: "*.Rds"
    publishDir "${params.outputDir}/signatures/${caller}/pdf", mode: 'copy', pattern: "*.pdf"

    input:
    set val(caller), val(type), val(sampleID), file(sample_vcf) from sample_sig_good_filtered

    output:
    file("${signatures_rds}")
    file("${signatures_plot_Rds}")
    file("${signatures_pieplot_Rds}")
    file("${signatures_plot_pdf}") into signatures_plots
    file("${signatures_pieplot_pdf}") into signatures_pie_plots
    set val("${caller}"), val("${type}"), file("${signatures_plot_pdf}") into signatures_plots_caller
    set val("${caller}"), val("${type}"), file("${signatures_pieplot_pdf}") into signatures_pieplots_caller
    set val(caller), val(type), val(sampleID), file("${signatures_weights_tsv}") into signatures_weights
    set val(sampleID), file("${signatures_rds}"), file("${signatures_plot_Rds}"), file("${signatures_pieplot_Rds}") into sample_signatures
    val(sampleID) into done_deconstructSigs_signatures

    script:
    prefix = "${sampleID}.${caller}.${type}"
    signatures_weights_tsv = "${prefix}.signatures.weights.tmp"
    signatures_rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_Rds']}"
    signatures_plot_pdf = "${prefix}.signatures.plot.pdf"
    signatures_plot_Rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_plot_Rds']}"
    signatures_pieplot_pdf = "${prefix}.signatures.pieplot.pdf"
    signatures_pieplot_Rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_pieplot_Rds']}"
    """
    deconstructSigs_make_signatures.R \
    "${sampleID}" \
    "${sample_vcf}" \
    "${signatures_rds}" \
    "${signatures_plot_pdf}" \
    "${signatures_plot_Rds}" \
    "${signatures_pieplot_pdf}" \
    "${signatures_pieplot_Rds}" \
    "${signatures_weights_tsv}"
    """
}

process update_signatures_weights {
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/signatures/${caller}/weights", mode: 'copy'

    input:
    set val(caller), val(type), val(sampleID), file(weights_tsv) from signatures_weights

    output:
    file("${output_tsv}") into updated_signatures_weights

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_tsv = "${prefix}.signatures.weights.tsv"
    """
    cat "${weights_tsv}" | \
    paste-col.py --header "Sample" -v "${sampleID}"  | \
    paste-col.py --header "VariantCallerType" -v "${type}"  | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${output_tsv}"
    """
}
updated_signatures_weights.collectFile(name: ".${signatures_weights_file}", keepHeader: true).set { updated_signatures_weights_collected }

process update_update_signatures_weights_collected {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from updated_signatures_weights_collected

    output:
    file("${output_file}") into all_signatures_weights

    script:
    output_file = "${signatures_weights_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}

// need to re-format the channel for combination with reporting channels
sample_signatures.map { sampleID, signatures_rds, signatures_plot_Rds, signatures_pieplot_Rds ->
    return([ sampleID, [ signatures_rds, signatures_plot_Rds, signatures_pieplot_Rds ] ])
}
.set { sample_signatures_reformated }

// group all plots by variant caller and type for merging
signatures_plots_caller.groupTuple(by: [0,1]).set{ signatures_plots_caller_grouped }
process merge_signatures_plots {
    // combine all signatures plots into a single PDF
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/signatures", mode: 'copy'

    input:
    set val(caller), val(type), file(input_files:'*') from signatures_plots_caller_grouped
    // file(input_files:'*') from signatures_plots.toList()

    output:
    file("${output_file}")
    val("${output_file}") into done_merge_signatures_plots

    when:
    input_files.size() > 0

    script:
    prefix = "${caller}.${type}"
    output_file="signatures.${prefix}.pdf"
    """
    gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile="${output_file}" ${input_files}
    """
}

// group all plots by variant caller and type for merging
signatures_pieplots_caller.groupTuple(by: [0,1]).set{ signatures_pieplots_caller_grouped }
process merge_signatures_pie_plots {
    // combine all signatures plots into a single PDF
    tag "${caller}.${type}"
    publishDir "${params.outputDir}/signatures", mode: 'copy'

    input:
    set val(caller), val(type), file(input_files:'*') from signatures_pieplots_caller_grouped

    output:
    file("${output_file}")
    val(output_file) into done_merge_signatures_pie_plots

    when:
    input_files.size() > 0

    script:
    prefix = "${caller}.${type}"
    output_file="signatures.pie.${prefix}.pdf"
    """
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${output_file}" ${input_files}
    """
}















// SETUP CHANNELS FOR PAIRED TUMOR-NORMAL STEPS FOR CNV analysis
// get all the combinations of samples & pairs
samples_bam3.combine(samples_pairs2) // [ sampleID, sampleBam, tumorID, normalID ]
                    .filter { item -> // only keep combinations where sample is same as tumor pair sample
                        def sampleID = item[0]
                        def sampleBam = item[1]
                        def tumorID = item[2]

                        sampleID == tumorID
                     }
                    .map { sampleID, sampleBam, tumorID, normalID -> // re-order the elements
                        def tumorBam = sampleBam
                        return [ tumorID, tumorBam, normalID ]
                    }
                    .combine(samples_bam4)  // combine again to get the samples & files again
                    .filter { item -> // keep only combinations where the normal ID matches the new sample ID
                        def tumorID = item[0]
                        def tumorBam = item[1]
                        def normalID = item[2]
                        def sampleID = item[3]
                        def sampleBam = item[4]
                        normalID == sampleID
                    }
                    .map {tumorID, tumorBam, normalID, sampleID, sampleBam  -> // re arrange the elements
                        def normalBam = sampleBam
                        def comparisonID = "${tumorID}_${normalID}"
                        return [ comparisonID, tumorID, tumorBam, normalID, normalBam ]
                    }//.subscribe { println "[samples_bam3]${it}"}
                    .tap { samples_pairs_bam_ch }
                    .combine(ref_fasta13) // add reference genome and targets
                    .combine(ref_fai13)
                    .combine(ref_dict13)
                    .combine(targets_annotated_bed)
                    .set {
                      samples_bam_ref_targets
                    }

process cnvkit {
    // CNV calling on the tumor-normal pairs
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), val(normalID), file(normalBam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_annotated_bed_file) from samples_bam_ref_targets

    output:
    file("${output_cns}")
    file("${output_finalcnr}")
    set val(comparisonID), val(tumorID), val(normalID), file("${output_cnr}"), file("${output_call_cns}"), file("${segment_gainloss}") into sample_cnvs

    script:
    prefix = "${comparisonID}"
    tumorBamID = "${tumorBam}".replaceFirst(/.bam$/, "")
    normal_reference = "normal_reference.cnn"
    tmp_cns = "${tumorBamID}.cns"
    tmp_cnr = "${tumorBamID}.cnr"
    output_cns = "${prefix}.cns"
    output_cnr = "${prefix}.cnr"
    output_finalcnr = "${prefix}.final.cnr"
    output_call_cns = "${prefix}.call.cns"
    segment_gainloss = "${prefix}.segment-gainloss.txt"
    """
    # running cnvkit pipeline on paried tumor/normal bams using batch mode, required tumorBam, normalBam, reference fasta and bed.
    cnvkit.py batch "${tumorBam}" \
    --normal "${normalBam}" \
    --targets "${targets_annotated_bed_file}" \
    --fasta "${ref_fasta}" \
    --output-reference "${normal_reference}" \
    -p \${NSLOTS:-\${NTHREADS:-1}}

    # Produces ${tmp_cns} and ${tmp_cnr}, rename to ${output_cns} and ${output_cnr}
    mv ${tmp_cns} ${output_cns}
    mv ${tmp_cnr} ${output_cnr}

    # Given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number.
    cnvkit.py call --filter cn "${output_cns}"

    # Identify targeted genes with copy number gain or loss above or below a threshold, -t :threshold and -m:number of bins
    cnvkit.py gainloss "${output_cnr}" -s "${output_call_cns}" -t 0.3 -m 5 > "${segment_gainloss}"
    cnvkit.py gainloss "${output_cnr}" -t 0.3 -m 5 > "${output_finalcnr}"

    """
}

// get a copy of the channel to use with the CNV pool
// combine against the CNV Pool file
samples_dd_bam5.combine(cnv_pool_ch)
.map { sampleID, bam, bai, cnv_pool_file ->
    // remap the channel to replace the Normal with CNV Pool
    def cnv_poolID = "CNV-Pool"
    def new_comparisonID = "${sampleID}_${cnv_poolID}"
    return([ new_comparisonID, sampleID, bam, cnv_poolID, cnv_pool_file ])
}
.set { samples_cnv_pool_ch }
// samples_cnv_pool_ch.subscribe { println "[samples_cnv_pool_ch] ${it}" }

process cnvkit_pooled_reference {
    // CNV calling on the tumor-poolednormal pairs
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(new_comparisonID), val(sampleID), file(bam), val(cnv_poolID), file(cnv_pool_file) from samples_cnv_pool_ch

    output:
    file("${output_cns}")
    file("${output_finalcnr}")
    set val(new_comparisonID), val(sampleID), val(cnv_poolID), file("${output_cnr}"), file("${output_call_cns}"), file("${segment_gainloss}") into sample_cnvs_pooledreference

    script:
    prefix = "${new_comparisonID}"
    tumorBamID = "${bam}".replaceFirst(/.bam$/, "")
    tmp_cns = "${tumorBamID}.cns"
    tmp_cnr = "${tumorBamID}.cnr"
    output_cns = "${prefix}.cns"
    output_cnr = "${prefix}.cnr"
    output_finalcnr = "${prefix}.final.cnr"
    output_call_cns = "${prefix}.call.cns"
    segment_gainloss = "${prefix}.segment-gainloss.txt"
    """
    # running cnvkit pipeline on tumor/pooled_normal reference using batch mode, required tumorBam, pooledReference.cnn.
    cnvkit.py batch "${bam}" \
    -r "${cnv_pool_file}" \
    -p \${NSLOTS:-\${NTHREADS:-1}}

    # Produces ${tmp_cns} and ${tmp_cnr}, rename to ${output_cns} and ${output_cnr}
    mv ${tmp_cns} ${output_cns}
    mv ${tmp_cnr} ${output_cnr}

    # Given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number.
    cnvkit.py call --filter cn "${output_cns}"

    # Identify targeted genes with copy number gain or loss above or below a threshold, -t :threshold and -m:number of bins
    cnvkit.py gainloss "${output_cnr}" -s "${output_call_cns}" -t 0.3 -m 5 > "${segment_gainloss}"
    cnvkit.py gainloss "${output_cnr}" -t 0.3 -m 5 > "${output_finalcnr}"
    """
}


// only keep files with at least 1 variant for TMB analysis
sample_cnvs_pooledreference.mix( sample_cnvs )
.filter { comparisonID, tumorID, normalID, cnr, call_cns, segment_gainloss ->
    def count = cnr.readLines().size()
    if (count <= 1) log.warn "${comparisonID} doesn't have enough lines in cnr and will not be processed"
    count > 1
}
.filter { comparisonID, tumorID, normalID, cnr, call_cns, segment_gainloss ->
    def count = segment_gainloss.readLines().size()
    if (count <= 1) log.warn "${comparisonID} doesn't have enough lines in segment_gainloss and will not be processed"
    count > 1
}
.set { sample_cnvs_filtered }

process cnvkit_gene_segments {
    // filter for segement and ratio genes that are in common between segment and ratio files; "trusted genes", avoid artifacts in the results
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(cnr), file(call_cns), file(segment_gainloss) from sample_cnvs_filtered

    output:
    set val(comparisonID), val(tumorID), val(normalID), file(segment_gainloss), file("${trusted_genes}") into sample_cnv_gene_segments

    script:
    prefix = "${comparisonID}"
    segment_genes = "${prefix}.segment-genes.txt"
    ratio_genes = "${prefix}.ratio-genes.txt"
    trusted_genes = "${prefix}.trusted-genes.txt"
    """
    # Generate segment genes (from .cns) and ratio genes (from .cnr) by applying threshold
    cnvkit.py gainloss "${cnr}" -s "${call_cns}" -t 0.3 -m 5 | tail -n+2 | cut -f1 | sort > "${segment_genes}"
    cnvkit.py gainloss "${cnr}" -t 0.3 -m 5 | tail -n+2 | cut -f1 | sort > "${ratio_genes}"

    # Take intersection of segment and ratio as trusted genes (final set of cnv genes)
    comm -12 "${ratio_genes}" "${segment_genes}" > "${trusted_genes}"
    """
}

// make sure that there are trusted genes present
sample_cnv_gene_segments.filter { comparisonID, tumorID, normalID, segment_gainloss, trusted_genes ->
    def count = trusted_genes.readLines().size()
    if (count <= 1) log.warn "${comparisonID} doesn't have enough lines in trusted_genes and will not be processed"
    count > 1
}.set { sample_cnv_gene_segments_filtered }

process cnvkit_extract_trusted_genes {
    // find gainloss segments with the trusted genes
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(segment_gainloss), file(trusted_genes) from sample_cnv_gene_segments_filtered

    output:
    set val(tumorID), val(normalID), file("${output_final_cns}") into cnvs_cns
    set val(comparisonID), val(tumorID), val(normalID), file("${output_final_cns}") into (cnvs_cns2, cnvs_cns3)


    script:
    prefix = "${comparisonID}"
    trusted_genes = "${prefix}.trusted-genes.txt"
    output_final_cns = "${prefix}.${filemap['cnvkit']['suffix']['final_cns']}" // use this file for report
    """
    # Get the trusted genes and their segment gain loss info from segment_gainloss file and write to final.cns
    cat "${segment_gainloss}" | head -n +1 > "${output_final_cns}"
    grep -w -f "${trusted_genes}" "${segment_gainloss}" >> "${output_final_cns}"
    """
}

process cnvkit_extract_trusted_genes_update {
    // convert flagstat output to a flat table
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(cns) from cnvs_cns2

    output:
    file("${output_file}") into cnvkit_extract_trusted_genes_updated


    script:
    prefix = "${comparisonID}"
    output_file = "${prefix}.${filemap['cnvkit']['suffix']['final_cns']}.sample.tsv"
    """
    paste-col.py -i "${cns}" --header "Tumor" -v "${tumorID}" | \
    paste-col.py --header "Normal" -v "${normalID}" | \
    paste-col.py --header "Comparison" -v "${comparisonID}" > \
    "${output_file}"
    """
}
cnvkit_extract_trusted_genes_updated.collectFile(name: ".cnv.tsv", keepHeader: true).set { cnvkit_extract_trusted_genes_collected }

process update_cnvkit_extract_trusted_genes_collected {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from cnvkit_extract_trusted_genes_collected

    output:
    file("${output_file}")

    script:
    output_file = "cnv.tsv"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
 }

process cnvkit_plotly {

    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(cns) from cnvs_cns3

    output:
    file("${output_html}")
    file("${output_pdf}")

    script:
    output_html = "${comparisonID}.cnvkit.plotly.html"
    output_pdf = "${comparisonID}.cnvkit.plotly.pdf"
    """
    cnvplot.R "${cns}" "${output_html}" "${output_pdf}"

    """
}

process snp_pileup {
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(snp_vcf), file(snp_vcf_tbi) from samples_dd_bam_noHapMap_pairs2.combine(common_snp_vcf).combine(common_snp_vcf_tbi)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file(output_cnvsnp) into snp_pileup

    script:
    caller = "FACETS"
    callerType = "cnv"
    prefix = "${comparisonID}.${caller}.${callerType}"
    output_cnvsnp = "${prefix}.snp_pileup.txt"
    """
    snp-pileup \
    -g \
    -q15 \
    -Q20 \
    -P100 \
    -r25,0 \
    "${snp_vcf}" \
    tmp.gz \
    "${normalBam}" \
    "${tumorBam}"

    # need to remove 'chr' from the table to resolve bugs in facets
    zcat tmp.gz | sed 's|chr||g' > "${output_cnvsnp}"
    """
}

// Need to filter out the samples that did not have enough entries for FACETS analysis
// need >1 entry (>2 lines in file)
snp_pileup_good = Channel.create()
snp_pileup_bad = Channel.create()
snp_pileup.choice( snp_pileup_good, snp_pileup_bad ){ items ->
    def caller = items[0]
    def callerType = items[1]
    def comparisonID = items[2]
    def tumorID = items[3]
    def normalID = items[4]
    def snp_pileup_txt = items[5]
    def output_ch = 1 // bad by default
    def num_lines = 0 // lines counter
    def required_lines = 2

    // make sure that the gz has at least 2 lines, then stop counting
    snp_pileup_txt.withReader { reader ->
            while (line = reader.readLine()) {
                if ( num_lines > required_lines ) {
                    output_ch = 0 // output to 'good' channel
                    break
                    }
                num_lines++
            }
        }
    return(output_ch)
}
snp_pileup_bad.map { caller, callerType, comparisonID, tumorID, normalID, snp_pileup_txt ->
    def reason = "Too few lines in snp_pileup_txt, skipping FACETS"
    def output = [comparisonID, caller, callerType, reason, "${snp_pileup_txt}"].join('\t')
    return(output)
}.set { snp_pileup_bad_logs }


process snp_pileup_check_variance {
    // need to check the variance of the given snp-pileup to make sure theres enough
    // data that FACETS will not break; needs to return some numeric value
    // bad data reurns 'NA'
    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file(snp_pileup_txt) from snp_pileup_good

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file(snp_pileup_txt), file("${output_file}") into snp_pileup_variances

    script:
    prefix = "${comparisonID}.${caller}.${callerType}"
    output_file = "${prefix}.variance.txt"
    """
    facets-check-variance.R "${snp_pileup_txt}" "${output_file}"
    """
}

// Need to filter out samples with a non-numeric variance value
snp_pileup_variance_good = Channel.create()
snp_pileup_variance_bad = Channel.create()
snp_pileup_variances.choice( snp_pileup_variance_good, snp_pileup_variance_bad ){ items ->
    def caller = items[0]
    def callerType = items[1]
    def comparisonID = items[2]
    def tumorID = items[3]
    def normalID = items[4]
    def snp_pileup_txt = items[5]
    def variance_txt = items[6]

    def output_ch = 1 // bad by default
    def line = new File("${variance_txt}").withReader { line = it.readLine() }

    if ( line != "NA" ){
        output_ch = 0
    }

    return(output_ch)
}
snp_pileup_variance_bad.map { caller, callerType, comparisonID, tumorID, normalID, snp_pileup_txt, variance_txt ->
    def reason = "SNP Pileup had invalid variance value, skipping FACETS"
    def output = [comparisonID, caller, callerType, reason, "${snp_pileup_txt}"].join('\t')
    return(output)
}.set { snp_pileup_variance_bad_logs }

process facets {
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file(snp_pileup_txt), file(snp_pileup_variance) from snp_pileup_variance_good

    output:
    file("${output_segment}")
    file("${output_pdf}")

    script:
    prefix = "${comparisonID}.${caller}.${callerType}"
    output_segment = "${prefix}.segment.csv"
    output_pdf = "${prefix}.plot.pdf"
    """
    facets.R "${snp_pileup_txt}" "${output_pdf}" "${output_segment}"
    """
}





// ~~~~~~~ QC: Homozygous SNP Comparison ~~~~~~ //
// compare the homozygous SNP's between sample pairs across unpaired variant calling
// overlap the two to make sure that the same are found in both

// [ caller, type, sampleID, annotations_tsv, tumorID, normalID ]
// combine against sample pairs to get tumor-normal pairs
annotations_annovar_tables3.combine(samples_pairs3)
    // filter for samples that match tumor ID
    .filter { items ->
        def caller = items[0]
        def type = items[1]
        def sampleID = items[2]
        def sample_tsv = items[3]
        def tumorID = items[4]
        def normalID = items[5]

        sampleID == tumorID
    }
    // rearrange
    .map { caller, type, sampleID, sample_tsv, tumorID, normalID ->
        def tumor_tsv = sample_tsv
        return [ caller, type, tumorID, tumor_tsv, normalID ]
    }
    // combine back against annotation tables to get normal annotation table
    .combine(annotations_annovar_tables4)
    // filter for only the ones that match the exact normal sample
    .filter { items ->
        def caller = items[0]
        def type = items[1]
        def tumorID = items[2]
        def tumor_tsv = items[3]
        def normalID = items[4]
        def sample_caller = items[5]
        def sample_type = items[6]
        def sampleID = items[7]
        def sample_tsv = items[8]

        def ID_match = normalID == sampleID
        def type_match = type == sample_type
        def caller_match = caller == sample_caller
        def all_matches = [ ID_match, type_match, caller_match ]
        for ( match in all_matches ){
            if ( match == false ){
                return false
            }
        }
        return true
    }
    .map { caller, type, tumorID, tumor_tsv, normalID, sample_caller, sample_type, sampleID, sample_tsv ->
        def normal_tsv = sample_tsv
        def comparisonID = "${tumorID}_${normalID}"
        return [ caller, type, comparisonID, tumorID, tumor_tsv, normalID, normal_tsv ]
    }
    .filter { caller, type, comparisonID, tumorID, tumor_tsv, normalID, normal_tsv ->
        // dont use indels
        type != "indel"
    }
    .set { annotations_tables_paired }
    // .subscribe { println "[annotations_annovar_tables4]: ${it}" }

process overlap_snp_filter {
    // need to pre-filter the annotation tables based on QC criteria
    publishDir "${params.outputDir}/snp_overlap/${caller}/filtered", mode: 'copy'

    input:
    set val(caller), val(type), val(comparisonID), val(tumorID), file(tumor_tsv), val(normalID), file(normal_tsv) from annotations_tables_paired

    output:
    set val(caller), val(type), val(comparisonID), val(tumorID), file("${tumor_output}"), val(normalID), file("${normal_output}") into annotations_tables_paired_filtered

    script:
    tumor_prefix = "${tumorID}.${caller}.${type}"
    tumor_output = "${tumor_prefix}.snp-overlap.filtered.tsv"
    normal_prefix = "${normalID}.${caller}.${type}"
    normal_output = "${normal_prefix}.snp-overlap.filtered.tsv"
    """
    snp-overlap-filter.py -c "${caller}" -i "${tumor_tsv}" | \
    paste-col.py --header "Sample" -v "${tumorID}" | \
    paste-col.py --header "Tumor" -v "${tumorID}" | \
    paste-col.py --header "Normal" -v "${normalID}" | \
    paste-col.py --header "VariantCallerType" -v "${type}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${tumor_output}"

    snp-overlap-filter.py -c "${caller}" -i "${normal_tsv}" | \
    paste-col.py --header "Sample" -v "${normalID}" | \
    paste-col.py --header "Tumor" -v "${tumorID}" | \
    paste-col.py --header "Normal" -v "${normalID}" | \
    paste-col.py --header "VariantCallerType" -v "${type}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${normal_output}"
    """
}

// Need to filter out the samples that did not have any remaining variants and log those entries
annotations_tables_paired_filtered_good = Channel.create()
annotations_tables_paired_filtered_bad = Channel.create()

// filter out samples with empty variant tables
annotations_tables_paired_filtered.choice( annotations_tables_paired_filtered_good, annotations_tables_paired_filtered_bad ){ items ->
    def caller = items[0]
    def type = items[1]
    def comparisonID = items[2]
    def tumorID = items[3]
    def tumor_tsv = items[4]
    def normalID = items[5]
    def normal_tsv = items[6]

    def output_ch = 1 // bad by default
    // long count = Files.lines(tsv).count()
    long tumor_count = Files.lines(tumor_tsv).count()
    long normal_count = Files.lines(normal_tsv).count()

    // good if either has >2 lines; 1 header and one entry
    if (tumor_count > 2) output_ch = 0
    if (normal_count > 2) output_ch = 0

    return(output_ch)
}

// log the bad sample combinations
annotations_tables_paired_filtered_bad.map { caller, type, comparisonID, tumorID, tumor_tsv, normalID, normal_tsv ->
    def reason = "No variants remaining after filtering, skipping homozygous SNP overlap"
    def output = [comparisonID, tumorID, normalID, type, caller, reason, "${tumor_tsv},${normal_tsv}"].join('\t')
    return(output)
}.set { annotations_tables_paired_filtered_bad_logs }

process overlap_snps {
    publishDir "${params.outputDir}/snp_overlap/${caller}/final", mode: 'copy'
    input:
    set val(caller), val(type), val(comparisonID), val(tumorID), file(tumor_tsv), val(normalID), file(normal_tsv) from annotations_tables_paired_filtered_good

    output:
    set val(caller), val(type), val(comparisonID), val(tumorID), val(normalID), file("${output_aggr_table}") into sample_snp_overlap_aggr

    script:
    prefix = "${comparisonID}.${caller}.${type}"
    output_matrix = "${prefix}.snp-overlap.matrix.tsv"
    output_table = "${prefix}.snp-overlap.tsv"
    output_aggr_table = "${prefix}.snp-overlap.aggregate.tsv"
    output_plot = "${prefix}.snp-overlap.pdf"
    """
    snp-overlap.R \
    "${tumor_tsv}" \
    "${normal_tsv}" \
    "${output_matrix}" \
    "${output_table}" \
    "${output_aggr_table}" \
    "${output_plot}"
    """
}

process update_overlap_snp_table {
    publishDir "${params.outputDir}/snp_overlap/${caller}", mode: 'copy'
    input:
    set val(caller), val(type), val(comparisonID), val(tumorID), val(normalID), file(tsv) from sample_snp_overlap_aggr

    output:
    file("${output_file}") into updated_snp_overlap_aggr

    script:
    prefix = "${comparisonID}.${caller}.${type}"
    output_file = "${prefix}.snp-overlap.aggregate.updated.tsv"
    """
    cat "${tsv}" | \
    paste-col.py --header "Tumor" -v "${tumorID}" | \
    paste-col.py --header "Normal" -v "${normalID}" | \
    paste-col.py --header "VariantCallerType" -v "${type}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > "${output_file}"
    """
}
updated_snp_overlap_aggr.collectFile(name: ".${snp_overlap_file}", keepHeader: true, storeDir: "${params.outputDir}").set { snp_overlap_collected }


process update_snp_overlap_collected {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from snp_overlap_collected

    output:
    file("${output_file}") into snp_overlap_collected_updated

    script:
    output_file = "${snp_overlap_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
 }




// ~~~~~~ QC: SeraCare selected variants comparison ~~~~~ //
// set val(caller), val(type), val(sampleID), file("${annotations_tsv}") into (annotations_annovar_tables2, annotations_annovar_tables3, annotations_annovar_tables4, annotations_annovar_tables5)
// seracare_sample_ids
annotations_annovar_tables5.combine(seracare_sample_ids).filter { items ->
    def caller = items[0]
    def type = items[1]
    def sampleID = items[2]
    def tsv = items[3]
    def seracare_sample_id = items[4]

    seracare_sample_id == sampleID
}
.filter { items ->
    // dont process 'indel' calling results
    def caller = items[0]
    def type = items[1]
    def sampleID = items[2]
    def tsv = items[3]
    def seracare_sample_id = items[4]

    ! "${type}".toLowerCase().contains("indel")
}
.filter { items ->
    // dont process 'indel' calling results
    def caller = items[0]
    def type = items[1]
    def sampleID = items[2]
    def tsv = items[3]
    def seracare_sample_id = items[4]

    caller == "LoFreq"
}
.map { caller, type, sampleID, tsv, seracare_sample_id ->
    return([ caller, type, sampleID, tsv ])
}
.combine(seracare_selected_tsv)
.set { seracare_annotations }

process filter_seracare_selected_variants {
    // filter for SeraCare Selected Variants and add extra annotation to the .tsv
    input:
    set val(caller), val(type), val(sampleID), file(tsv), file(selected_tsv) from seracare_annotations

    output:
    set val(caller), val(type), val(sampleID), file("${output_file}") into seracare_variants

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_file = "${prefix}.seracare-selected.tsv"
    """
    seracare-selected-variant-filter.py -i "${tsv}" -o "${output_file}" --seracare "${selected_tsv}"
    """
}

process confidence_intervals_seracare {
    // need to calculate the confidence interval for variant calling for the mutations from SeraCare samples
    input:
    set val(caller), val(type), val(sampleID), file(tsv) from seracare_variants

    output:
    file("${output_file}") into seracare_variants_ci

    script:
    prefix = "${sampleID}.${caller}.${type}"
    output_file = "${prefix}.seracare-selected.ci.tsv"
    """
    variant-confidence-intervals.R "${tsv}" "${output_file}" "${SeraCareErrorRate}"
    """
}

process collect_seracare_annotation_tables {
    // combine all variants into a single table

    input:
    file('t*') from seracare_variants_ci.collect()

    output:
    file("${output_file}") into collected_seracare_variants

    script:
    output_file = "all_seracare_annotations.tmp.tsv"
    """
    concat-tables.py t* > "${output_file}"
    """
}

process update_collect_seracare_annotation_tables {
    // add labels to the table to output
    publishDir "${params.outputDir}/annotations", mode: 'copy'

    input:
    file(table) from collected_seracare_variants

    output:
    file("${output_file}") into collected_updated_seracare_variants

    script:
    output_file = "${all_seracare_annotations_file}" // annotations.SeraCare.tsv
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" | \
    paste-col.py --header "GitBranch" -v "${params.GIT_CURRENT_BRANCH}" | \
    paste-col.py --header "GitTag" -v "${params.GIT_CURRENT_TAG}" > \
    "${output_file}"
    """
}
collected_updated_seracare_variants.mix(seracare_selected_tsv2).set { seracare_pool_items } // to report





// ~~~~~ IGV Snapshot ~~~~~ //
// visualization of select variants in IGV

vcf_tsv_pairs2.filter { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, tsv ->
    // only keep entries with >1 line (header)
    def count = tsv.readLines().size()
    count > 1
}.set { vcf_tsv_pairs_filtered }

process igv_filter_variant {
    // filter the tsv variant table to find variants to snapshot

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(tsv) from vcf_tsv_pairs_filtered

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${output_file}") into igv_filtered_variants

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    output_file = "${prefix}.igv-annotations.tsv"
    """
    igv-variant-filter.py -c "${caller}" -i "${tsv}" > "${output_file}"
    """
}

igv_filtered_variants.filter { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, tsv ->
    // only keep entries with >1 line (header)
    def count = tsv.readLines().size()
    count > 1
}.set { igv_refiltered_variants }

process igv_tsv_to_bed {
    // convert the variants tsv table into a bed file
    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(tsv) from igv_refiltered_variants

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${output_file}") into igv_regions

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    output_file = "${prefix}.igv-regions.bed"
    """
    variant-tsv2bed.py -i "${tsv}" -o "${output_file}"
    """
}

// TODO: Need to chunk the bed file because some of them are too big to run, need to parallelize more
process igv_bed_line_chunk {
    // split the .bed file into new files based on desired lines in each file

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(bed) from igv_regions

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file('*') into igv_regions_split

    script:
    numLines = 30 // IGV snaphot takes ~2s per target region; 30 regions = ~1min execution
    """
    split-bed-lines.py "${bed}" "${numTargetSplitLines}"
    """
}

// 'bed_files' may be a single file or a list of files;
// need to convert all to list, then emit all files individually, and get the chunk label
igv_regions_split.map { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed_files ->
    // check if its a list or not
    def is_list = bed_files instanceof Collection
    if(is_list){
        return([ caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed_files ])
    } else {
        // make it a list
        def bed_list = [ bed_files ]
        return([ caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed_list ])
    }
}
.transpose() // emit each bed chunk individually
.map { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed ->
    // get number at the end of the file basename to denote the chunk Label
    def bedChunkLabel = "${bed.name}".findAll(/\d*$/)[0]

    return([ caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed, bedChunkLabel ])
}
.set { igv_regions_chunked }

// Combine against the tumor-normal pairs bam files
igv_regions_chunked.combine(samples_dd_ra_rc_bam_noHapMap_pairs)
.filter { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed, bedChunkLabel, paired_comparisonID, paired_tumorID, tumorBam, tumorBai, paired_normalID, normalBam, normalBai ->
    def comparisonID_match = comparisonID ==  paired_comparisonID
    def tumorID_match = tumorID == paired_tumorID
    def normalID_math = normalID == paired_normalID
    comparisonID_match && tumorID_match && normalID_math
}.map { caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed, bedChunkLabel, paired_comparisonID, paired_tumorID, tumorBam, tumorBai, paired_normalID, normalBam, normalBai ->
    return([ caller, callerType, comparisonID, tumorID, normalID, chunkLabel, bed, bedChunkLabel, tumorBam, tumorBai, normalBam, normalBai ])
}
.set { igv_regions_bams }

process igv_snapshot {
    // run IGV headlessly to make snapshots

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(bed), val(bedChunkLabel), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from igv_regions_bams

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), val(bedChunkLabel), file("${snapshotDirectory}") into igv_snaphots

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}.${bedChunkLabel}"
    batchscript = "snapshots.bat"
    snapshotDirectory = "${prefix}"
    image_height = 750
    genome = 'hg19'
    """
    mkdir "${snapshotDirectory}"

    # make IGV snapshot batch script
    make-igv-batchscript.py \
    "${tumorBam}" "${normalBam}" \
    -r "${bed}" \
    -d "${snapshotDirectory}" \
    -b "${batchscript}" \
    --height "${image_height}" \
    --genome "${genome}"

    # run IGV headlessly with xvfb
    xvfb-run \
    --auto-servernum \
    --server-num=1 \
    igv.sh \
    -b "${batchscript}"
    """
}

// group all the snapshot dirs together by category
igv_snaphots.groupTuple(by: [0, 1, 2, 3, 4, 5])
.set { igv_snaphots_grouped }

process aggregate_snapshots_bedChunkLabels {
    stageInMode "copy"
    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), val(bedChunkLabels), file("*") from igv_snaphots_grouped

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${snapshotDirectory}") into aggregated_snapshots

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    snapshotDirectory = "${prefix}"
    """
    mkdir "${snapshotDirectory}"
    find . -type f -name "*.png" -exec mv {} "${snapshotDirectory}/" \\;
    """
}

aggregated_snapshots.groupTuple(by: [0, 1, 2, 3, 4])
.set { reaggregated_snapshots }

process aggregate_snapshots_chunkLabels {
    stageInMode "copy"
    publishDir "${params.outputDir}/igv-snapshots", mode: 'copy'
    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabels), file("*") from reaggregated_snapshots

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), file("${snapshotDirectory}")

    script:
    prefix = "${comparisonID}.${caller}.${callerType}"
    snapshotDirectory = "${prefix}"
    """
    mkdir "${snapshotDirectory}"
    find . -type f -name "*.png" -exec mv {} "${snapshotDirectory}/" \\;
    """
}

// ~~~~~~~~ REPORTING ~~~~~~~ //
// collect from all processes to make sure they are finished before starting reports
done_copy_samplesheet.concat(
    done_print_metadata,
    done_fastq_merge,
    done_trimmomatic,
    done_fastqc_trim,
    done_alignment,
    done_sambamba_flagstat,
    done_sambamba_dedup,
    done_sambamba_dedup_flagstat,
    done_qc_target_reads_gatk_genome,
    done_qc_target_reads_gatk_pad500,
    done_qc_target_reads_gatk_pad100,
    done_qc_target_reads_gatk_bed,
    done_gatk_RealignerTargetCreator,
    done_gatk_IndelRealigner,
    done_gatk_BaseRecalibrator,
    done_gatk_BaseRecalibratorBQSR,
    done_gatk_PrintReads,
    done_lofreq,
    done_gatk_hc,
    done_delly2,
    done_deconstructSigs_signatures,
    done_merge_signatures_plots,
    done_merge_signatures_pie_plots,
    done_msisensor,
    done_mutect2,
    done_annotate,
    done_annotate_pairs,
    done_collect_annotation_tables,
    done_sambamba_dedup_log_table,
    done_sambamba_flagstat_table,
    done_sambamba_dedup_flagstat_table,
    done_update_coverage_tables,
    done_update_interval_tables,
    done_eval_sample_vcf,
    done_eval_pair_vcf,
    done_samtools_flagstat_table_update,
    done_update_sambamba_dedup_log_table,
    done_update_samtools_dedup_flagstat_table,
    done_update_updated_coverage_tables_collected,
    done_update_collect_annotation_tables,
    done_update_updated_coverage_interval_tables_collected
    )
    .set { all_done }


// collect failed log messages
failed_samples.mix(samples_vcfs_tsvs_bad_logs, sample_sig_bad_logs)
    .collectFile(name: "failed.tsv", storeDir: "${params.outputDir}", newLine: true)
    .set { failed_log_ch }
failed_pairs.mix(pairs_vcfs_tsvs_bad_logs, annotations_tables_paired_filtered_bad_logs, snp_pileup_bad_logs, snp_pileup_variance_bad_logs)
    .collectFile(name: "failed.pairs.tsv", storeDir: "${params.outputDir}", newLine: true)
    .set{ failed_pairs_log_ch }


// need a channel for items that may or may not have been produced depending on the conditional nature of the pipeline config
// e.g. paired vs unpaired runs, some files might not exist
Channel.fromPath('.placeholder1')
.mix(snp_overlap_collected_updated, seracare_pool_items)
.set { report_inputs }

process custom_analysis_report {
    // create a batch report for all samples in the analysis
    publishDir "${params.outputDir}", mode: 'copy'
    stageInMode "copy"

    input:
    file(report_items: '*') from analysis_report_files.collect()
    file(all_annotations_file) from all_annotations_file_ch
    file(samplesheet_output_file) from samplesheet_output_file_ch
    file(sample_coverage_file) from sample_coverage_file_ch
    file(interval_coverage_file) from interval_coverage_file_ch
    file(meta_file) from metadata_ch
    file(reads_dedup_table) from sambamba_dedup_log_table_ch
    file(flagstat_table) from samtools_flagstat_table_ch
    file(dedup_flagstat_table) from samtools_dedup_flagstat_table_ch
    file(targets_metrics_table) from targets_metrics_ch
    file(failed_log) from failed_log_ch
    file(failed_pairs_log) from failed_pairs_log_ch
    file(targets_annotations_file) from annotated_targets
    file(signatures_weights) from all_signatures_weights
    file(tmb_file) from tmbs_collected
    file(git_json) from git_json_ch
    file(extra_report_inputs: "*") from report_inputs.collect()

    output:
    file("${html_output}")

    script:
    prefix = "${runID}"
    html_output = "${prefix}.report.html"
    """
    # convert report file symlinks to copies of original files, because knitr doesnt work well unless all report files are in pwd
    for item in *.Rmd *.css *.bib; do
        if [ -L "\${item}" ]; then
            sourcepath="\$(python -c "import os; print(os.path.realpath('\${item}'))")"
            echo ">>> resolving source file: \${sourcepath}"
            rsync -va "\${sourcepath}" "\${item}"
        fi
    done

    R --vanilla <<E0F
    rmarkdown::render(input = "main.Rmd",
    params = list(
        annotations_file = "${all_annotations_file}",
        samplesheet_file = "${samplesheet_output_file}",
        sample_coverage_file = "${sample_coverage_file}",
        interval_coverage_file = "${interval_coverage_file}",
        meta_file = "${meta_file}",
        reads_dedup_table = "${reads_dedup_table}",
        flagstat_table = "${flagstat_table}",
        dedup_flagstat_table = "${dedup_flagstat_table}",
        targets_metrics_table = "${targets_metrics_table}",
        failed_log = "${failed_log}",
        failed_pairs_log = "${failed_pairs_log}",
        targets_annotations_file = "${targets_annotations_file}",
        signatures_weights_file = "${signatures_weights}",
        git_json_file = "${git_json}",
        snp_overlap_file = "${snp_overlap_file}",
        seracare_annotations_file = "${all_seracare_annotations_file}",
        seracare_selected_variants = "${SeraCareSelectedTsvFile}"
        ),
    output_format = "html_document",
    output_file = "${html_output}")
    E0F
    """
}

// channel for sample reports; gather items per-sample from other processes
sampleIDs.map { sampleID ->
    // dummy file to pass through channel
    def placeholder = file(".placeholder1")

    return([ sampleID, placeholder ])
}
// add items from other channels (unpaired steps)
.concat(sample_signatures_reformated) // [ sampleID, [ sig_file1, sig_file2, ... ] ]
// group all the items by the sampleID, first element in each
.groupTuple() // [ sampleID, [ [ sig_file1, sig_file2, ... ], .placeholder, ... ] ]
// need to flatten any nested lists
.map { sampleID, fileList ->
    def newFileList = fileList.flatten()

    return([ sampleID, newFileList ])
}
.into { sample_output_files; sample_output_files2 }

// channel for items from sample pairs
samples_pairs_with_NA.map { tumorID, normalID ->
    // dummy file to pass through channel
    def placeholder = file(".placeholder2")

    return([ tumorID, normalID, placeholder ])
}
// add tumor-normal process output channels here
.concat(cnvs_cns)
// group by tumor and normal id's
.groupTuple(by: [0,1])
.map { tumorID, normalID, fileList ->
    // un-nest any nested lists of files
    def newFileList = fileList.flatten()

    return([ tumorID, normalID, newFileList ])
}
.set { sample_pairs_output_files }

// combine the channels
sample_pairs_output_files.combine(sample_output_files2)
.filter { tumorID, normalID, tumorNormalFileList, sampleID, sampleFileList ->
    sampleID == tumorID
}
.map { tumorID, normalID, tumorNormalFileList, sampleID, sampleFileList ->
    return([ tumorID, normalID, tumorNormalFileList, sampleFileList ])
}
.set { sample_pairs_output_files2 }
// .subscribe { println "[sample_output_files2] ${it}" }

process custom_sample_report {
    // create per-sample reports
    publishDir "${params.outputDir}/reports", mode: 'copy'

    input:
    set val(tumorID), val(normalID), file(tumorNormalFiles: "*"), file(sampleFiles: "*"), file(report_items: '*') from sample_pairs_output_files2.combine(samples_report_files)

    output:
    file("${html_output}")

    script:
    prefix = "${tumorID}"
    html_output = "${prefix}.report.html"
    """
    # convert report file symlinks to copies of original files, because knitr doesnt work well unless all report files are in pwd
    for item in *.Rmd *.css *.bib; do
        if [ -L "\${item}" ]; then
            sourcepath="\$(python -c "import os; print(os.path.realpath('\${item}'))")"
            echo ">>> resolving source file: \${sourcepath}"
            rsync -va "\${sourcepath}" "\${item}"
        fi
    done

    R --vanilla <<E0F
    rmarkdown::render(input = "main.Rmd",
    params = list(
        sampleID = "${tumorID}",
        signatures_Rds = "${filemap['deconstructSigs']['suffix']['signatures_Rds']}",
        signatures_plot_Rds = "${filemap['deconstructSigs']['suffix']['signatures_plot_Rds']}",
        signatures_pieplot_Rds = "${filemap['deconstructSigs']['suffix']['signatures_pieplot_Rds']}",
        paired_normal = "${normalID}"
        ),
    output_format = "html_document",
    output_file = "${html_output}")
    E0F
    """
}

Channel.fromPath("${params.outputDir}").set { output_dir_ch }
process multiqc {
    // automatic reporting based on detected output; might take a while to parse and create report
    publishDir "${params.outputDir}/qc", mode: 'copy'

    input:
    val(all_vals) from all_done.collect()
    file(output_dir) from output_dir_ch

    output:
    file("${output_html}")
    file("${output_data}")

    when:
    disable_multiqc != true

    script:
    prefix = "${runID}"
    multiqc_html = "multiqc_report.html"
    multiqc_data = "multiqc_data"
    output_html = "${prefix}.multiqc.html"
    output_data = "multiqc-data"
    """
    multiqc "${output_dir}"
    mv "${multiqc_html}" "${output_html}"
    mv "${multiqc_data}" "${output_data}"
    """
}


// ~~~~~~~~~~~~~~~ PIPELINE COMPLETION EVENTS ~~~~~~~~~~~~~~~~~~~ //
// gather email file attachments
// Channel.fromPath( file(samplesheet) ).set{ samples_analysis_sheet }
// def attachments = samples_analysis_sheet.toList().getVal()

workflow.onComplete {

    def status = "NA"

    if(workflow.success) {
        status = "SUCCESS"

    } else {
        status = "FAILED"
    }

    def log_message = "${status} [${workflowTimestamp}] [${workflow.complete.format('dd-MMM-yyyy HH:mm:ss')}] [${workflow.duration}] [${workflow.runName}] [${workflow.sessionId}] [${workflow.scriptId ?: '-'}] [${workflow.nextflow.version} ${workflow.nextflow.build}]\n"

    nxf_completion_log = new File(params.nxf_completion_log)
    nxf_completion_log.append(log_message)

    def msg = """
        Pipeline execution summary
        ---------------------------
        Success           : ${workflow.success}
        exit status       : ${workflow.exitStatus}
        Launch time       : ${workflow.start.format('dd-MMM-yyyy HH:mm:ss')}
        Ending time       : ${workflow.complete.format('dd-MMM-yyyy HH:mm:ss')} (duration: ${workflow.duration})
        Launch directory  : ${workflow.launchDir}
        Work directory    : ${workflow.workDir.toUriString()}
        Project directory : ${workflow.projectDir}
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow container: ${workflow.container ?: '-'}
        container engine  : ${workflow.containerEngine?:'-'}
        Nextflow run name : ${workflow.runName}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})


        The command used to launch the workflow was as follows:

        ${workflow.commandLine}

        --
        This email was sent by Nextflow
        cite doi:10.1038/nbt.3820
        http://nextflow.io
        """
        .stripIndent()

    sendMail {
        to "${params.email_to}"
        from "${params.email_from}"
        // attach attachments // attachments keep breaking hold off on these
        subject "[${params.workflowLabel}] ${status}"

        body
        """
        ${msg}
        """
        .stripIndent()
    }
}
