import java.nio.file.Files;
import java.text.SimpleDateFormat;
import groovy.json.JsonSlurper;

//
// NGS580 Target exome analysis for 580 gene panel
//

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

// ~~~~~~~~~~ CONFIGURATION ~~~~~~~~~~ //
// configure pipeline settings
// overriden by nextflow.config and CLI args
params.configFile = "config.json"
params.outputDir = "output"
params.reportDir = "report"
params.targetsBed = "targets.bed"
params.targetsrefFlatBed = "targets.refFlat.580.bed"
params.probesBed = "probes.bed"
params.samplesheet = null
params.runID = null
def username = params.username // from nextflow.config
def workflowTimestamp = "${workflow.start.format('yyyy-MM-dd-HH-mm-ss')}"
def outputDirPath = new File(params.outputDir).getCanonicalPath()
def reportDirPath = new File(params.reportDir).getCanonicalPath()
String localhostname = java.net.InetAddress.getLocalHost().getHostName();
// Date now = new Date()
// SimpleDateFormat timestamp = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss")
// currentDirPath = new File(System.getProperty("user.dir")).getCanonicalPath()

// load the JSON config, if present
def jsonSlurper = new JsonSlurper()
def PipelineConfig
def PipelineConfigFile_obj = new File("${params.configFile}")
if ( PipelineConfigFile_obj.exists() ) {
    log.info("Loading configs from ${params.configFile}")
    String PipelineConfigJSON = PipelineConfigFile_obj.text
    PipelineConfig = jsonSlurper.parseText(PipelineConfigJSON)
}

// check for Run ID
// 0. use CLI passed arg
// 1. check for config.json values
// 2. use the name of the current directory
def projectDir = "${workflow.projectDir}"
def projectDirname = new File("${projectDir}").getName()
def runID
if( params.runID == null ){
    if ( PipelineConfig && PipelineConfig.containsKey("runID") && PipelineConfig.runID != null ) {
        log.info("Loading runID from ${params.configFile}")
        runID = "${PipelineConfig.runID}"
    } else {
        runID = "${projectDir}"
    }
} else {
    runID = "${projectDir}"
}

// Check for samplesheet;
// 0. Use CLI passed samplesheet
// 1. check for config.json values
// 2. Check for samples.analysis.tsv in current directory
def default_samplesheet = "samples.analysis.tsv"
def default_samplesheet_obj = new File("${default_samplesheet}")
def samplesheet
if(params.samplesheet == null){
    if ( PipelineConfig && PipelineConfig.containsKey("samplesheet") && PipelineConfig.samplesheet != null ) {
        log.info("Loading samplesheet from ${params.configFile}")
        samplesheet = PipelineConfig.samplesheet
    } else if( default_samplesheet_obj.exists() ){
        samplesheet = default_samplesheet_obj.getCanonicalPath()
    } else {
        log.error("No samplesheet found, please provide one with '--samplesheet'")
        exit 1
    }
} else {
    samplesheet = params.samplesheet
}

// Enable or disable some pipeline steps here TODO: better config management for this
disable_multiqc = true // for faster testing of the rest of the pipeline
disable_msisensor = true // breaks on very small demo datasets
disable_delly2 = true

// names of some important output files to use throughout the pipeline
def all_annotations_file = "annotations.tsv"
def samplesheet_output_file = 'samplesheet.tsv'
def sample_coverage_file = "coverage.samples.tsv"
def interval_coverage_file = "coverage.intervals.tsv"
def signatures_weights_file = "signatures.weights.tsv"

// load a mapping dict to use for keeping track of the names and suffixes for some files throughout the pipeline
String filemapJSON = new File("filemap.json").text
def filemap = jsonSlurper.parseText(filemapJSON)

// ~~~~~ START WORKFLOW ~~~~~ //
log.info "~~~~~~~ NGS580 Pipeline ~~~~~~~"
log.info "* Launch time:        ${workflowTimestamp}"
log.info "* Run ID:             ${runID}"
log.info "* Samplesheet:        ${samplesheet}"
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
Channel.fromPath( file(params.targetsBed) ).set{ targets_bed }
Channel.fromPath( file(params.targetsrefFlatBed) ).set{ targets_refFlat_bed }

// reference files
Channel.fromPath( file(params.targetsBed) ).into { targets_bed; targets_bed2; targets_bed3; targets_bed4; targets_bed5; targets_bed6; targets_bed7; targets_bed8; targets_bed9 }

Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2; ref_fasta3; ref_fasta4; ref_fasta5; ref_fasta6; ref_fasta7; ref_fasta8; ref_fasta9; ref_fasta10; ref_fasta11; ref_fasta12; ref_fasta13 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2; ref_fai3; ref_fai4; ref_fai5; ref_fai6; ref_fai7; ref_fai8; ref_fai9; ref_fai10; ref_fai11; ref_fai12; ref_fai13 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2; ref_dict3; ref_dict4; ref_dict5; ref_dict6; ref_dict7; ref_dict8; ref_dict9; ref_dict10; ref_dict11; ref_dict12; ref_dict13 }

Channel.fromPath( file(params.ref_chrom_sizes) ).set{ ref_chrom_sizes }
Channel.fromPath( file(params.trimmomatic_contaminant_fa) ).set{ trimmomatic_contaminant_fa }
Channel.fromPath( file(params.ref_fa_bwa_dir) ).set{ ref_fa_bwa_dir }

Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf) ).into{ gatk_1000G_phase1_indels_vcf; gatk_1000G_phase1_indels_vcf2; gatk_1000G_phase1_indels_vcf3; gatk_1000G_phase1_indels_vcf4 }
Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf_idx) ).into{ gatk_1000G_phase1_indels_vcf_idx; gatk_1000G_phase1_indels_vcf_idx2; gatk_1000G_phase1_indels_vcf_idx3; gatk_1000G_phase1_indels_vcf_idx4 }

Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf) ).into{ mills_and_1000G_gold_standard_indels_vcf; mills_and_1000G_gold_standard_indels_vcf2; mills_and_1000G_gold_standard_indels_vcf3; mills_and_1000G_gold_standard_indels_vcf4 }
Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf_idx) ).into{ mills_and_1000G_gold_standard_indels_vcf_idx; mills_and_1000G_gold_standard_indels_vcf_idx2; mills_and_1000G_gold_standard_indels_vcf_idx3; mills_and_1000G_gold_standard_indels_vcf_idx4 }

Channel.fromPath( file(params.dbsnp_ref_vcf) ).into{ dbsnp_ref_vcf; dbsnp_ref_vcf2; dbsnp_ref_vcf3; dbsnp_ref_vcf4; dbsnp_ref_vcf5; dbsnp_ref_vcf6; dbsnp_ref_vcf7; dbsnp_ref_vcf8 }
Channel.fromPath( file(params.dbsnp_ref_vcf_idx) ).into{ dbsnp_ref_vcf_idx; dbsnp_ref_vcf_idx2; dbsnp_ref_vcf_idx3; dbsnp_ref_vcf_idx4; dbsnp_ref_vcf_idx5; dbsnp_ref_vcf_idx6; dbsnp_ref_vcf_idx7; dbsnp_ref_vcf_idx8 }

Channel.fromPath( file(params.cosmic_ref_vcf) ).into{ cosmic_ref_vcf; cosmic_ref_vcf2 }
Channel.fromPath( file(params.cosmic_ref_vcf_idx) ).into{ cosmic_ref_vcf_idx; cosmic_ref_vcf_idx2 }

Channel.fromPath( file(params.microsatellites) ).set{ microsatellites }
Channel.fromPath( file(params.ANNOVAR_DB_DIR) ).into { annovar_db_dir; annovar_db_dir2; annovar_db_dir3 }

// report and output dir
Channel.fromPath("${outputDirPath}").into { analysis_output; analysis_output2 }

// load analysis report files
Channel.fromPath("${reportDirPath}/analysis/*")
        .set { analysis_report_files }

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
        .into { samples_pairs; samples_pairs2 }

// read sample IDs from analysis sheet
Channel.fromPath( file(samplesheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sampleID = row['Sample']
            return(sampleID)
        }
        .unique()
        .into { sampleIDs; sampleIDs2 }

Channel.fromPath( file(samplesheet) ).set { samples_analysis_sheet }

// logging channels
Channel.from("Sample\tProgram\tNote\tFiles").set { failed_samples }
Channel.from("Comparison\tTumor\tNormal\tChrom\tProgram\tNote\tFiles").set { failed_pairs }

// ~~~~~ PIPELINE TASKS ~~~~~ //
// PREPROCESSING
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

process fastq_merge {
    // merge multiple R1 and R2 fastq files (e.g. split by lane) into a single fastq each
    tag "${sampleID}"

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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    paste-col.py --header "System" -v "${localhostname}" > \
    "${output_file}"
    """
}

process sambamba_dedup {
    // deduplicate alignments
    tag "${sampleID}"
    publishDir "${params.outputDir}/alignments/deduplicated", mode: 'copy'

    input:
    set val(sampleID), file(sample_bam) from samples_bam2

    output:
    set val(sampleID), file("${bam_file}") into samples_dd_bam, samples_dd_bam2
    file("${bai_file}")
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
    tag "${sampleID}"
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
    paste-col.py --header "System" -v "${localhostname}" > \
    "${output_file}"
    """
}

process samtools_dedup_flagstat {
    // dedup alignment stats
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    paste-col.py --header "System" -v "${localhostname}" > \
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
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
    tag "${sampleID}"
    publishDir "${params.outputDir}/alignments/stats", mode: 'copy'

    input:
    set val(sampleID), file(table1), file(table2), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict) from recalibrated_bases_table2_comb

    output:
    file("${csv_file}")
    file("${pdf_file}")
    val(sampleID) into done_gatk_AnalyzeCovariates

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
    tag "${sampleID}"
    publishDir "${params.outputDir}/alignments/recalibrated", mode: 'copy'

    input:
    set val(sampleID), file(table1), file(ra_bam_file), file(ra_bai_file), file(ref_fasta), file(ref_fai), file(ref_dict) from recalibrated_bases_table1_ref

    output:
    set val(sampleID), file("${ra_rc_bam_file}"), file("${ra_rc_bai_file}") into samples_dd_ra_rc_bam, samples_dd_ra_rc_bam2, samples_dd_ra_rc_bam3, samples_dd_ra_rc_bam4
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
                            samples_dd_ra_rc_bam_ref8 }


samples_dd_ra_rc_bam_ref.combine( dbsnp_ref_vcf3 )
                        .combine(dbsnp_ref_vcf_idx3)
                        .into { samples_dd_ra_rc_bam_ref_dbsnp;
                                samples_dd_ra_rc_bam_ref_dbsnp2 }



process qc_target_reads_gatk_genome {
    // calculate metrics on coverage across whole genome
    tag "${prefix}"
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
    -ct 10 -ct 50 -ct 100 -ct 500 \
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
    tag "${prefix}"
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
    -ct 10 -ct 50 -ct 100 -ct 500 \
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
    tag "${prefix}"
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
    -ct 10 -ct 50 -ct 100 -ct 500 \
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
    tag "${prefix}"
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
    -ct 10 -ct 50 -ct 100 -ct 500 \
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
    tag "${prefix}"
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
    paste-col.py --header "System" -v "${localhostname}" > \
    "${output_file}"
    """
}

process coverage_intervals_to_table {
    // convert the intervals coverage file into a formatted table for downstream usage
    tag "${prefix}"

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
    tag "${prefix}"

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
    --onetranscript \
    --outfile "${prefix}"

    merge-interval-tables.R "${interval_table}" "${annovar_output_txt}" "${avinput_file}" "${annotations_tsv}"
    """
}

process update_interval_tables {
    // add more information to the intervals coverage output
    tag "${prefix}"
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
    paste-col.py --header "System" -v "${localhostname}" > \
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
    tag "${sampleID}"
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_bgz_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from samples_dd_ra_rc_bam_ref_dbsnp

    output:
    set val(caller), val(sampleID), file("${filtered_vcf}"), file("${reformat_tsv}") into samples_lofreq_vcf
    set val(caller), val(sampleID), file("${filtered_vcf}") into samples_lofreq_vcf2
    file("${vcf_file}")
    file("${norm_vcf}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    file("${tsv_file}")
    file("${reformat_tsv}")
    file("${filtered_vcf}")
    val(sampleID) into done_lofreq

    script:
    caller = "LoFreq"
    prefix = "${sampleID}.${caller}"
    vcf_file = "${prefix}.vcf"
    vcf_bgz_file = "${prefix}.vcf.bgz"
    norm_vcf = "${prefix}.norm.vcf"
    filtered_vcf = "${prefix}.filtered.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    """
    lofreq call-parallel \
    --call-indels \
    --pp-threads \${NSLOTS:-\${NTHREADS:-1}} \
    --ref "${ref_fasta}" \
    --bed "${targets_bed_file}" \
    --out "${vcf_file}" \
    "${sample_bam}"

    bgzip -c "${vcf_file}" > "${vcf_bgz_file}"

    bcftools index "${vcf_bgz_file}"

    cat ${vcf_file} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"

    # do not report if frequency is less than 1%
    gatk.sh -T SelectVariants \
    -R "${ref_fasta}" \
    -V "${norm_vcf}" \
    -select "AF > 0.01"  \
    > "${filtered_vcf}"

    gatk.sh -T VariantsToTable \
    -R "${ref_fasta}" \
    -V "${filtered_vcf}" \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F SB -F INDEL -F CONSVAR -F HRUN \
    -o "${tsv_file}"

    # reformat and adjust the TSV table for consistency downstream
    # add extra columns to the VCF TSV file for downstream
    reformat-vcf-table.py -c LoFreq -s "${sampleID}" -i "${tsv_file}" | \
    paste-col.py --header "Sample" -v "${sampleID}"  | \
    paste-col.py --header "VariantCaller" -v "${caller}" > \
    "${reformat_tsv}"
    """
}

process gatk_hc {
    // variant calling
    tag "${sampleID}"
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"


    input:
    set val(sampleID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx) from samples_dd_ra_rc_bam_ref_dbsnp2

    output:
    set val(caller), val(sampleID), file("${filtered_vcf}") into sample_vcf_hc
    set val(caller), val(sampleID), file("${filtered_vcf}"), file("${reformat_tsv}") into sample_vcf_hc2
    set val(caller), val(sampleID), file("${filtered_vcf}") into sample_vcf_hc3
    file("${vcf_file}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    file("${norm_vcf}")
    file("${tsv_file}")
    file("${reformat_tsv}")
    file("${filtered_vcf}")
    val(sampleID) into done_gatk_hc

    script:
    caller = "HaplotypeCaller"
    prefix = "${sampleID}.${caller}"
    vcf_file = "${prefix}.vcf"
    norm_vcf = "${prefix}.norm.vcf"
    filtered_vcf = "${prefix}.filtered.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
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

    # report if
    # alternate allele freq (allele depth / depth) greater than 0.05 ; 5%
    # more than 5 variant call supporting reads
    # quality reads present (reported depth >0)
    gatk.sh -T SelectVariants \
    -R "${ref_fasta}" \
    -V "${norm_vcf}" \
    --sample_name "${sampleID}" \
    -select "vc.getGenotype('${sampleID}').getAD().1 / vc.getGenotype('${sampleID}').getDP() > 0.05" \
    -select "vc.getGenotype('${sampleID}').getAD().1 > 5" \
    -select "vc.getGenotype('${sampleID}').getDP() > 0" \
    > "${filtered_vcf}"

    gatk.sh -T VariantsToTable \
    -R "${ref_fasta}" \
    -V "${filtered_vcf}" \
    -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN \
    -GF AD -GF DP \
    -o "${tsv_file}"

    # reformat and adjust the TSV table for consistency downstream
    # add extra columns to the VCF TSV file for downstream
    reformat-vcf-table.py -c HaplotypeCaller -s "${sampleID}" -i "${tsv_file}" | \
    paste-col.py --header "Sample" -v "${sampleID}" | \
    paste-col.py --header "VariantCaller" -v "${caller}" > \
    "${reformat_tsv}"
    """
}

sample_vcf_hc3.mix(samples_lofreq_vcf2)
                .combine(dbsnp_ref_vcf4)
                .combine(dbsnp_ref_vcf_idx4)
                .combine(ref_fasta4)
                .combine(ref_fai4)
                .combine(ref_dict4)
                .set { samples_filtered_vcfs }
process eval_sample_vcf {
    // calcaulte variant metrics
    tag "${sampleID}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${eval_file}"

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(ref_fasta), file(ref_fai), file(ref_dict4) from samples_filtered_vcfs

    output:
    file("${eval_file}")
    val("${sampleID}") into done_eval_sample_vcf

    script:
    prefix = "${sampleID}.${caller}"
    eval_file = "${prefix}.eval.grp"
    """
    gatk.sh -T VariantEval \
    -R "${ref_fasta}" \
    -o "${eval_file}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --eval "${sample_vcf}"
    """
}


//
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
    tag "${sampleID}-${type}"
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


// Genomic Signatures
sample_vcf_hc_good = Channel.create()
sample_vcf_hc_bad = Channel.create()
deconstructSigs_variant_min = 55

sample_vcf_hc.choice( sample_vcf_hc_good, sample_vcf_hc_bad ){ items ->
    // make sure there are enough variants in the VCF to proceed!
    def caller = items[0]
    def sampleID = items[1]
    def filtered_vcf = items[2]
    def line_count = 0
    def num_variants = 0
    def output = 1 // bad by default
    def enough_variants = false // bad by default

    // count number of variants in the vcf file
    filtered_vcf.withReader { reader ->
        while (line = reader.readLine()) {
            if (!line.startsWith("#")) num_variants++
            if (num_variants > deconstructSigs_variant_min) {
                enough_variants = true
                break
                }
            line_count++
        }
    }

    if ( enough_variants==false ) {
        output = 1
    } else if ( enough_variants==true ) {
        output = 0
    }
    return(output)
}

sample_vcf_hc_bad.map {  caller, sampleID, filtered_vcf ->
    def reason = "Fewer than ${deconstructSigs_variant_min} variants in .vcf file, skipping genomic signatures"
    def output = [sampleID, caller, reason, filtered_vcf].join('\t')
    return(output)
}.set { sample_vcf_hc_bad_logs }

process deconstructSigs_signatures {
    // search for mutation signatures
    tag "${sampleID}"
    publishDir "${params.outputDir}/signatures/${caller}/Rds", mode: 'copy', pattern: "*.Rds"
    publishDir "${params.outputDir}/signatures/${caller}/pdf", mode: 'copy', pattern: "*.pdf"

    input:
    set val(caller), val(sampleID), file(sample_vcf) from sample_vcf_hc_good

    output:
    file("${signatures_rds}")
    file("${signatures_plot_Rds}")
    file("${signatures_pieplot_Rds}")
    file("${signatures_plot_pdf}") into signatures_plots
    file("${signatures_pieplot_pdf}") into signatures_pie_plots
    set val(caller), val(sampleID), file("${signatures_weights_tsv}") into signatures_weights
    set val(sampleID), file("${signatures_rds}"), file("${signatures_plot_Rds}"), file("${signatures_pieplot_Rds}") into sample_signatures
    val(sampleID) into done_deconstructSigs_signatures

    script:
    prefix = "${sampleID}.${caller}"
    signatures_weights_tsv = "${prefix}.signatures.weights.tmp"
    signatures_rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_Rds']}"
    signatures_plot_pdf = "${prefix}.signatures.plot.pdf"
    signatures_plot_Rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_plot_Rds']}"
    signatures_pieplot_pdf = "${prefix}.signatures.pieplot.pdf"
    signatures_pieplot_Rds = "${prefix}.${filemap['deconstructSigs']['suffix']['signatures_pieplot_Rds']}"
    """
    deconstructSigs_make_signatures.R "${sampleID}" "${sample_vcf}" "${signatures_rds}" "${signatures_plot_pdf}" "${signatures_plot_Rds}" "${signatures_pieplot_pdf}" "${signatures_pieplot_Rds}" "${signatures_weights_tsv}"
    """
}

process update_signatures_weights {
    tag "${sampleID}"
    publishDir "${params.outputDir}/signatures/${caller}/weights", mode: 'copy'

    input:
    set val(caller), val(sampleID), file(weights_tsv) from signatures_weights

    output:
    file("${output_tsv}") into updated_signatures_weights

    script:
    prefix = "${sampleID}.${caller}"
    output_tsv = "${prefix}.signatures.weights.tsv"
    """
    cat "${weights_tsv}" | \
    paste-col.py --header "Sample" -v "${sampleID}"  | \
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
    file("${output_file}")

    script:
    output_file = "${signatures_weights_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" > \
    "${output_file}"
    """
}
// need to re-format the channel for combination with reporting channels
sample_signatures.map { sampleID, signatures_rds, signatures_plot_Rds, signatures_pieplot_Rds ->
    return([ sampleID, [ signatures_rds, signatures_plot_Rds, signatures_pieplot_Rds ] ])
}
.set { sample_signatures_reformated }

process merge_signatures_plots {
    // combine all signatures plots into a single PDF
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(input_files:'*') from signatures_plots.toList()

    output:
    file("${output_file}")
    val("${output_file}") into done_merge_signatures_plots

    when:
    input_files.size() > 0

    script:
    output_file="signatures.pdf"
    """
    gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile="${output_file}" ${input_files}
    """
}


process merge_signatures_pie_plots {
    // combine all signatures plots into a single PDF
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(input_files:'*') from signatures_pie_plots.toList()

    output:
    file("${output_file}")
    val(output_file) into done_merge_signatures_pie_plots

    when:
    input_files.size() > 0

    script:
    output_file="signatures.pie.pdf"
    """
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${output_file}" ${input_files}
    """
}


// // REQUIRES ANNOTATIONS FOR DBSNP FILTERING
// process vaf_distribution_plot {
//     tag { "${sampleID}" }
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
// samples_dd_ra_rc_bam2 // set val(sampleID), file("${sampleID}.dd.ra.rc.bam"), file("${sampleID}.dd.ra.rc.bam.bai")
// samples_pairs // [ tumorID, normalID ]
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
                    .tap { samples_dd_ra_rc_bam_pairs } // make a channel for just this set of data
                    .combine(ref_fasta3) // add reference genome and targets
                    .combine(ref_fai3)
                    .combine(ref_dict3)
                    .combine(targets_bed4)
                    .tap {  samples_dd_ra_rc_bam_pairs_ref; // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) ]
                            samples_dd_ra_rc_bam_pairs2;
                            samples_dd_ra_rc_bam_pairs_ref3;
                            samples_dd_ra_rc_bam_pairs_ref4 }
                    .combine(microsatellites) // add MSI ref
                    .set { samples_dd_ra_rc_bam_pairs_ref_msi }




// get the unique chromosomes in the targets bed file
//  for per-chrom paired variant calling
Channel.fromPath( params.targetsBed )
            .splitCsv(sep: '\t')
            .map{row ->
                row[0]
            }
            .unique()
            .set{ chroms }

// add the reference .vcf
samples_dd_ra_rc_bam_pairs_ref.combine(dbsnp_ref_vcf2)
                            .combine(dbsnp_ref_vcf_idx2)
                            .combine(cosmic_ref_vcf2)
                            .combine(cosmic_ref_vcf_idx2)
                            .tap { samples_dd_ra_rc_bam_pairs_ref_gatk } // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), dbsnp, cosmic ]
                            // add the chroms
                            .combine( chroms ) // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), dbsnp, cosmic, chrom ]
                            .set { samples_dd_ra_rc_bam_pairs_ref_gatk_chrom }


// REQUIRES PAIRED SAMPLES BAM FILES
process msisensor {
    // detection of microsatellite instability regions
    // disable this since it always breaks on small sample datasets
    // this program is buggy;
    // TODO: find a method to pre-filter based on the number of alignments ??
    tag "${prefix}"
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
    tag "${prefix}"
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy', pattern: "*${vcf_file}"
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy', pattern: "*${norm_vcf}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${multiallelics_stats}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${realign_stats}"


    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(cosmic_ref_vcf), file(cosmic_ref_vcf_idx), val(chrom) from samples_dd_ra_rc_bam_pairs_ref_gatk_chrom

    output:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file("${norm_vcf}") into vcfs_mutect2
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file("${norm_vcf}") into samples_mutect3 // to VariantEval eval_pair_vcf
    file("${vcf_file}")
    file("${norm_vcf}")
    file("${multiallelics_stats}")
    file("${realign_stats}")
    val(comparisonID) into mutect2_sampleIDs
    val(comparisonID) into done_mutect2

    script:
    caller = "MuTect2"
    prefix = "${comparisonID}.${chrom}.${caller}"
    bed_subset = "${prefix}.bed"
    vcf_file = "${prefix}.vcf"
    norm_vcf = "${prefix}.norm.vcf"
    filtered_vcf = "${prefix}.filtered.vcf"
    multiallelics_stats = "${prefix}.bcftools.multiallelics.stats.txt"
    realign_stats = "${prefix}.bcftools.realign.stats.txt"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    """
    # subset the target regions for the given chromosome
    subset_bed.py "${chrom}" "${targets_bed}" > "${bed_subset}"

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
    --intervals "${bed_subset}" \
    --interval_padding 10 \
    --input_file:tumor "${tumorBam}" \
    --input_file:normal "${normalBam}" \
    --out "${vcf_file}"

    # normalize and split vcf entries
    cat ${vcf_file} | \
    bcftools norm --multiallelics -both --output-type v - 2>"${multiallelics_stats}" | \
    bcftools norm --fasta-ref "${ref_fasta}" --output-type v - 2>"${realign_stats}" > \
    "${norm_vcf}"
    """
}

process filter_vcf_pairs {
    // filter the .vcf for tumor-normal pairs
    tag "${prefix}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"

    input:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from vcfs_mutect2.combine(ref_fasta11).combine(ref_fai11).combine(ref_dict11)

    output:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file("${filtered_vcf}") into filtered_vcf_pairs // to tsv

    script:
    prefix = "${comparisonID}.${chrom}.${caller}"
    filtered_vcf = "${prefix}.filtered.vcf"
    if( caller == 'MuTect2' )
        """
        # filter VCF
        # report if:
        # only keep 'PASS' entries

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :

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
        """
    else
        error "Invalid caller: ${caller}"
}

process vcf_to_tsv_pairs {
    tag "${prefix}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcf_pairs.combine(ref_fasta12).combine(ref_fai12).combine(ref_dict12)

    output:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file(vcf), file("${reformat_tsv}") into vcf_tsv_pairs // to annotation

    script:
    prefix = "${comparisonID}.${chrom}.${caller}"
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

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}"  | \
        paste-col.py --header "Tumor" -v "${tumorID}"  | \
        paste-col.py --header "Normal" -v "${normalID}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"
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
    tag "${prefix}"
    publishDir "${params.outputDir}/variants/${caller}/stats", mode: 'copy', pattern: "*${eval_file}"

    input:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file(pairs_vcf), file(dbsnp_ref_vcf), file(dbsnp_ref_vcf_idx), file(ref_fasta), file(ref_fai), file(ref_dict4) from pairs_filtered_vcfs

    output:
    file("${eval_file}")
    val("${comparisonID}") into done_eval_pair_vcf

    script:
    prefix = "${comparisonID}.${chrom}.${caller}"
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

samples_lofreq_vcf.concat(sample_vcf_hc2).set { samples_vcfs_tsvs }

// filter out samples with empty variant tables
samples_vcfs_tsvs.choice( samples_vcfs_tsvs_good, samples_vcfs_tsvs_bad ){ items ->
    def caller =  items[0]
    def sampleID =  items[1]
    def sample_vcf =  items[2]
    def sample_tsv =  items[3]
    def output = 1 // bad by default
    long count = Files.lines(sample_tsv).count()
    if (count > 1) output = 0 // good if has >1 lines
    return(output)
}

samples_vcfs_tsvs_bad.map { caller, sampleID, sample_vcf, sample_tsv ->
    def reason = "Too few lines in sample_tsv, skipping annotation"
    def output = [sampleID, caller, reason, "${sample_vcf},${sample_tsv}"].join('\t')
    return(output)
}.set { samples_vcfs_tsvs_bad_logs }

vcf_tsv_pairs.choice( pairs_vcfs_tsvs_good, pairs_vcfs_tsvs_bad ){ items ->
    def caller = items[0]
    def comparisonID = items[1]
    def tumorID = items[2]
    def normalID = items[3]
    def chrom = items[4]
    def sample_vcf = items[5]
    def sample_tsv = items[6]
    def output = 1 // bad by default
    long count = Files.lines(sample_tsv).count()
    if (count > 1) output = 0 // good if has >1 lines
    return(output)
}

pairs_vcfs_tsvs_bad.map { caller, comparisonID, tumorID, normalID, chrom, sample_vcf, sample_tsv ->
    def reason = "Too few lines in sample_tsv, skipping annotation"
    def output = [comparisonID, tumorID, normalID, chrom, caller, reason, "${sample_vcf},${sample_tsv}"].join('\t')
    return(output)
}.set { pairs_vcfs_tsvs_bad_logs }

process annotate {
    // annotate variants
    tag "${prefix}"
    publishDir "${params.outputDir}/annotations/${caller}", mode: 'copy', pattern: "*${annotations_tsv}"

    input:
    set val(caller), val(sampleID), file(sample_vcf), file(sample_tsv), file(annovar_db_dir) from samples_vcfs_tsvs_good.combine(annovar_db_dir)

    output:
    file("${annotations_tsv}") into annotations_tables
    set val(sampleID), val(caller), file("${annotations_tsv}") into annotations_annovar_tables
    file("${avinput_file}")
    file("${avinput_tsv}")
    val(sampleID) into done_annotate

    script:
    prefix = "${sampleID}.${caller}"
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
        --onetranscript \
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
        --onetranscript \
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
    tag "${prefix}"
    publishDir "${params.outputDir}/annotations/${caller}", mode: 'copy', pattern: "*${annotations_tsv}"

    input:
    set val(caller), val(comparisonID), val(tumorID), val(normalID), val(chrom), file(sample_vcf), file(sample_tsv), file(annovar_db_dir) from pairs_vcfs_tsvs_good.combine(annovar_db_dir2)

    output:
    file("${annotations_tsv}") into annotations_tables_pairs
    file("${avinput_file}")
    file("${avinput_tsv}")
    file("${annovar_output_txt}")
    file("${annovar_output_vcf}")
    val(comparisonID) into done_annotate_pairs

    script:
    prefix = "${comparisonID}.${chrom}.${caller}"
    avinput_file = "${prefix}.avinput"
    avinput_tsv = "${prefix}.avinput.tsv"
    annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
    annovar_output_vcf = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.vcf"
    annotations_tsv = "${prefix}.annotations.tsv"
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
        --onetranscript \
        --outfile "${prefix}"

        # get values from .avinput file
        printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${avinput_tsv}"
        cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"

        # merge the tables together
        merge-vcf-tables.R "${sample_tsv}" "${annovar_output_txt}" "${avinput_tsv}" "${annotations_tsv}"
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

process update_collect_annotation_tables {
    // add labels to the table to output
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file(table) from collected_annotation_tables

    output:
    file("${output_file}") into all_annotations_file_ch
    val('') into done_update_collect_annotation_tables

    script:
    output_file = "${all_annotations_file}"
    """
    paste-col.py -i "${table}" --header "Run" -v "${runID}" | \
    paste-col.py --header "Time" -v "${workflowTimestamp}" | \
    paste-col.py --header "Session" -v "${workflow.sessionId}" | \
    paste-col.py --header "Workflow" -v "${workflow.runName}" | \
    paste-col.py --header "Location" -v "${workflow.projectDir}" | \
    paste-col.py --header "System" -v "${localhostname}" > \
    "${output_file}"
    """
}


// ~~~~ Tumor Burden Analysis ~~~~ //
samples_dd_ra_rc_bam4.combine(ref_fasta10)
.combine(ref_fai10)
.combine(ref_dict10)
.combine(targets_bed8)
.set { samples_bam_ref }

process gatk_CallableLoci {
    tag "${prefix}"
    publishDir "${params.outputDir}/loci", mode: 'copy'

    input:
    set val(sampleID), file(bamfile), file(baifile), file(genomeFa), file(genomeFai), file(genomeDict), file(targetsBed) from samples_bam_ref

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
    --intervals "${targetsBed}" \
    -o "${output_bed}"

    grep -E 'CALLABLE|PASS' "${output_bed}" > "${output_bed_pass}" || touch "${output_bed_pass}" # exit code 1 if no matches
    """
}

process callable_loci_table {
    tag "${prefix}"
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
annotations_annovar_tables.filter { sampleID, caller, anno_tsv ->
    def count = anno_tsv.readLines().size()
    count > 1
}.set { annotations_annovar_tables_filtered }

process tmb_filter_variants {
    tag "${prefix}"
    publishDir "${params.outputDir}/tmb/${caller}/Rdata", mode: 'copy', pattern: "*${output_Rdata}"
    publishDir "${params.outputDir}/tmb/${caller}/variants", mode: 'copy', pattern: "*${output_variants}"

    input:
    set val(sampleID), val(caller), file(anno_tsv) from annotations_annovar_tables_filtered

    output:
    file("${output_Rdata}")
    file("${output_variants}") into tmb_filtered_variants
    set val(sampleID), val(caller), file("${output_variants}") into tmb_filtered_variants2

    script:
    prefix = "${sampleID}.${caller}"
    output_Rdata = "${prefix}.tmb_filter_variants.Rdata"
    output_variants = "${prefix}.annotations.tmb_filtered.tsv"
    tmpFile = "tmp"
    """
    tmb-variant-filter.R "${output_Rdata}" "${tmpFile}" "${anno_tsv}"
    cat "${tmpFile}" | \
    paste-col.py --header "Sample" -v "${sampleID}" | \
    paste-col.py --header "Caller" -v "${caller}" > "${output_variants}"
    """
}

loci_tables2.map{ sampleID, summary_table, callable_txt ->
    // read the number of callable loci from the file, just pass the value
    def callable_loci = callable_txt.readLines()[0]
    return([sampleID, summary_table, "${callable_loci}"])
}
// combine against the filtered tmb variants
.combine(tmb_filtered_variants2)
.filter{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, filtered_annotations ->
    // only keep the matches between channels
    sampleID_anno == sampleID_loci
}
.map{ sampleID_loci, loci_summary_table, callable_loci, sampleID_anno, caller, filtered_annotations ->
    // read the number of variants from the file
    def num_variants = filtered_annotations.readLines().size() - 1
    return([sampleID_loci, caller, callable_loci, num_variants])
}
.filter{ sampleID_loci, caller, callable_loci, num_variants ->
    // only keep samples that had callable loci
    callable_loci > 0
}
.set { loci_annotations }

process calculate_tmb {
    tag "${prefix}"
    publishDir "${params.outputDir}/tmb/${caller}/values", mode: 'copy', pattern: "*${output_tmb}"

    input:
    set val(sampleID), val(caller), val(loci), val(variants) from loci_annotations

    output:
    file("${output_tmb}") into tmbs

    script:
    prefix = "${sampleID}.${caller}"
    output_tmb = "${prefix}.tmb.tsv"
    """
    tmb=\$( calc-tmb.py ${variants} ${loci} )
    printf 'SampleID\tCaller\tnBases\tnVariants\tTMB\n' > "${output_tmb}"
    printf "${sampleID}\t${caller}\t${loci}\t${variants}\t\${tmb}\n" >> "${output_tmb}"
    """
}
tmbs.collectFile(name: 'tmb.tsv', keepHeader: true, storeDir: "${params.outputDir}")


// SETUP CHANNELS FOR PAIRED TUMOR-NORMAL STEPS FOR CNV analysis
// samples_pairs // [ tumorID, normalID ]
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
                    .combine(ref_fasta13) // add reference genome and targets
                    .combine(ref_fai13)
                    .combine(ref_dict13)
                    .combine(targets_refFlat_bed)
                    .set {
                      samples_bam_ref_targets
                    }

process cnvkit {
    publishDir "${params.outputDir}/cnv", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), val(normalID), file(normalBam), file(ref_fasta13), file(ref_fai13), file(ref_dict13), file(targets_refFlat_bed) from samples_bam_ref_targets

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
    --targets "${targets_refFlat_bed}" \
    --fasta "${ref_fasta13}" \
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

// only keep files with at least 1 variant for TMB analysis
sample_cnvs.filter { comparisonID, tumorID, normalID, cnr, call_cns, segment_gainloss ->
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
    input:
    set val(comparisonID), val(tumorID), val(normalID), file(segment_gainloss), file(trusted_genes) from sample_cnv_gene_segments_filtered

    output:
    set val(tumorID), val(normalID), file("${output_final_cns}") into cnvs_cns

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
    done_gatk_AnalyzeCovariates,
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
failed_samples.concat(samples_vcfs_tsvs_bad_logs, sample_vcf_hc_bad_logs)
    .collectFile(name: "failed.tsv", storeDir: "${params.outputDir}", newLine: true)
    .set { failed_log_ch }
failed_pairs.concat(pairs_vcfs_tsvs_bad_logs)
    .collectFile(name: "failed.pairs.tsv", storeDir: "${params.outputDir}", newLine: true)
    .set{ failed_pairs_log_ch }


process custom_analysis_report {
    // create a batch report for all samples in the analysis
    tag "${html_output}"
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
        failed_pairs_log = "${failed_pairs_log}"
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
    tag "${prefix}"
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

    if(params.pipeline_email) {
        sendMail {
            to "${params.email_to}"
            from "${params.email_from}"
            // attach attachments // attachments keep breaking hold off on these
            subject "[${params.workflow_label}] ${status}"

            body
            """
            ${msg}
            """
            .stripIndent()
        }
    }
}
