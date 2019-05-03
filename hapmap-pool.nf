params.samplesheet = "samples.hapmap.tsv"
params.outputDir = "output-hapmap-pool"
params.outputBamName = "HapMap-pool"

def samplesheet = params.samplesheet
def workflowTimestamp = "${workflow.start.format('yyyy-MM-dd-HH-mm-ss')}"
def outputDirPath = new File(params.outputDir).getCanonicalPath()

// ~~~~~ START WORKFLOW ~~~~~ //
log.info "~~~~~~~ NGS580 HapMap Pool Pipeline ~~~~~~~"
log.info "* Launch time:        ${workflowTimestamp}"
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


// read samples from analysis samplesheet
Channel.fromPath( file(samplesheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sampleID = row['Sample']
            def bam = row['Bam']
            return([ sampleID, file(bam) ])
        }.set { samplesIDs_bams }

process fix_header {
    stageInMode "copy"
    input:
    set val(sampleID), file(bam) from samplesIDs_bams

    output:
    file("${output_file}") into fixed_header_bams

    script:
    output_file = "${sampleID}.fixed.bam"
    """
    # need to fix the headers in each .bam file; change the sample name to "HapMap-pool"
    # https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    samtools view -H "${bam}" | \
    sed 's|\\(^@RG.*\\)\\(SM:[^\\s]*\\)\\(.*\$\\)|\\1SM:${params.outputBamName}\\3|' | \
    samtools reheader - "${bam}" > "${output_file}"
    """
}

process bam_merge {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    file("*") from fixed_header_bams.collect()

    output:
    file("${output_bam}")
    file("${output_bai}")

    script:
    output_bam = "${params.outputBamName}.bam"
    output_bai = "${output_bam}.bai"
    """
    samtools merge --threads=\${NSLOTS:-\${NTHREADS:-1}} - *.bam | \
    samtools sort --threads=\${NSLOTS:-\${NTHREADS:-1}} > "${output_bam}"
    samtools index "${output_bam}"
    """
}
