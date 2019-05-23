params.samplesheet = "samples.cnv.tsv"
params.outputDir = "output-cnv-pool"
params.targetsBed = "targets/targets.annotated.580.bed"
params.targetsLabel = "580"
params.outputLabel = "CNV-Pool"
params.outputFileName = "${params.outputLabel}.${params.targetsLabel}"

def samplesheet = params.samplesheet
def workflowTimestamp = "${workflow.start.format('yyyy-MM-dd-HH-mm-ss')}"
def outputDirPath = new File(params.outputDir).getCanonicalPath()

// ~~~~~ START WORKFLOW ~~~~~ //
log.info "~~~~~~~ NGS580 CNV Pool Pipeline ~~~~~~~"
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
            def bam = row['Bam']
            return([ file(bam) ])
        }.set { sample_bams }

// reference files
Channel.fromPath( file("${params.targetsBed}") ).set { targets_bed }
Channel.fromPath( file(params.ref_fa) ).into { ref_fasta_ch; ref_fasta_ch2 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai_ch; ref_fai_ch2 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict_ch; ref_dict_ch2 }


sample_bams.combine(targets_bed)
    .combine(ref_fasta_ch)
    .combine(ref_fai_ch)
    .combine(ref_dict_ch)
    .set { sample_cnv_inputs }

process cnv_reference {
    // https://cnvkit.readthedocs.io/en/v0.9.0/pipeline.html
    // Make target and antitarget coverage.cnn for each normal sample
    input:
    set file(bam), file(bed), file(ref_fasta), file(ref_fai), file(ref_dict) from sample_cnv_inputs

    output:
    file("${output_file}") into cnns_ch

    script:
    output_file = "${bam}.cnn"
    """
    cnvkit.py batch \
    -n "${bam}" \
    --output-reference "${output_file}" \
    --targets "${bed}" \
    --fasta "${ref_fasta}" \
    -p \${NSLOTS:-\${NTHREADS:-1}}
    """
}
cnns_ch.toList().toList()
.combine(ref_fasta_ch2)
.combine(ref_fai_ch2)
.combine(ref_dict_ch2)
.set { cnn_ref_ch }

process cnv_pooledreference {
    // Combine all .cnn as PooledReference.cnn
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    set file("*.cnn"), file(ref_fasta), file(ref_fai), file(ref_dict) from cnn_ref_ch

    output:
    file("${output_file}")

    script:
    output_file = "${params.outputFileName}.cnn"
    """
    cnvkit.py reference *.cnn -f "${ref_fasta}" -o "${output_file}"
    """
}
