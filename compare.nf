params.old_unpaired_annotations = null
params.new_unpaired_annotations = null
params.old_paired_annotations = null
params.new_paired_annotations = null
params.outputDir = "output"
params.reportDir = "report"
def old_unpaired_annotations = new File("${params.old_unpaired_annotations}")
def new_unpaired_annotations = new File("${params.new_unpaired_annotations}")
def old_paired_annotations = new File("${params.old_paired_annotations}")
def new_paired_annotations = new File("${params.new_paired_annotations}")
def reportDirPath = new File(params.reportDir).getCanonicalPath()
def outputDirPath = new File(params.outputDir).getCanonicalPath()
def workflowTimestamp = "${workflow.start.format('yyyy-MM-dd-HH-mm-ss')}"

if( ! old_unpaired_annotations.exists() | ! new_unpaired_annotations.exists() | ! old_paired_annotations.exists() | ! new_paired_annotations.exists() | "${params.old_unpaired_annotations}" == "${params.new_unpaired_annotations}" | "${params.old_paired_annotations}" == "${params.new_paired_annotations}" ){
    log.error "Invalid files; old_unpaired_annotations: ${params.old_unpaired_annotations} , new_unpaired_annotations: ${params.new_unpaired_annotations}, new_paired_annotations: ${new_paired_annotations} , old_paired_annotations: ${old_paired_annotations}"
    exit 1
}

log.info "~~~~~~~ NGS580 Annotation Comparison ~~~~~~~"
log.info "* Launch time:        ${workflowTimestamp}"
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

Channel.fromPath("${reportDirPath}/compare/*")
        .into { report_files; report_files2 }
Channel.fromPath("${old_unpaired_annotations}").set{ old_unpaired_annotations_ch }
Channel.fromPath("${new_unpaired_annotations}").set{ new_unpaired_annotations_ch }
Channel.fromPath("${old_paired_annotations}").set{ old_paired_annotations_ch }
Channel.fromPath("${new_paired_annotations}").set{ new_paired_annotations_ch }

Channel.from(new File("${old_unpaired_annotations}").getCanonicalPath()).set{ old_unpaired_annotations_path_ch }
Channel.from(new File("${new_unpaired_annotations}").getCanonicalPath()).set{ new_unpaired_annotations_path_ch }
Channel.from(new File("${old_paired_annotations}").getCanonicalPath()).set{ old_paired_annotations_path_ch }
Channel.from(new File("${new_paired_annotations}").getCanonicalPath()).set{ new_paired_annotations_path_ch }


process compare_report {
    publishDir "${outputDirPath}", mode: 'copy'
    stageInMode "copy"

    input:
    set file("old.unpaired.annot.tsv"),
        file("new.unpaired.annot.tsv"),
        file("old.paired.annot.tsv"),
        file("new.paired.annot.tsv"),
        val(old_unpaired_annot_path),
        val(new_unpaired_annot_path),
        val(old_paired_annot_path),
        val(new_paired_annot_path) from old_unpaired_annotations_ch.combine(new_unpaired_annotations_ch)
                                                .combine(old_paired_annotations_ch)
                                                .combine(new_paired_annotations_ch)
                                                .combine(old_unpaired_annotations_path_ch)
                                                .combine(new_unpaired_annotations_path_ch)
                                                .combine(old_paired_annotations_path_ch)
                                                .combine(new_paired_annotations_path_ch)
    file(report_items: '*') from report_files.collect()

    output:
    file("${html_output}")
    file("old.unpaired.annot.tsv")
    file("new.unpaired.annot.tsv")
    file("old.paired.annot.tsv")
    file("new.paired.annot.tsv")
    file("*.tsv")

    script:
    html_output = "comparison.html"
    """
    R --vanilla <<E0F
    rmarkdown::render(
        input = "compare.Rmd",
    params = list(
        old_unpaired_annot = "old.unpaired.annot.tsv",
        new_unpaired_annot = "new.unpaired.annot.tsv",
        old_paired_annot = "old.paired.annot.tsv",
        new_paired_annot = "new.paired.annot.tsv",
        old_unpaired_annot_path = "${old_unpaired_annot_path}",
        new_unpaired_annot_path = "${new_unpaired_annot_path}",
        old_paired_annot_path = "${old_paired_annot_path}",
        new_paired_annot_path = "${new_paired_annot_path}"
    ),
    output_format = "html_document",
    output_file = "${html_output}"
    )
    E0F
    """
}
