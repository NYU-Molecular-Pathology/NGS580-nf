// NGS580 Target exome analysis for 580 gene panel

// pipeline settings; overriden by nextflow.config and CLI args
params.output_dir = "output"

// summary collectFile's
params.qc_coverage_gatk_file_basename = "qc_coverage_gatk.csv"
params.annotations_mutect2_file_basename = "annotations-mutect2.txt"
params.annotations_insertions_Delly2_file_basename = "annotations-insertions-Delly2.txt"
params.annotations_translocations_Delly2_file_basename = "annotations-translocations-Delly2.txt"
params.annotations_inversion_Delly2_file_basename = "annotations-inversions-Delly2.txt"
params.annotations_duplications_Delly2_file_basename = "annotations-duplications-Delly2.txt"
params.annotations_deletions_Delly2_file_basename = "annotations-deletions-Delly2.txt"
params.annotations_hc_file_basename = "annotations-hc.txt"
params.annotations_lofreq_file_basename = "annotations-lofreq.txt"


//
// DATA INPUT CHANNELS
//
// targets .bed file
Channel.fromPath( file(params.targets_bed) ).set{ targets_bed }

// reference files
Channel.fromPath( file(params.targets_bed) ).into { targets_bed; targets_bed2; targets_bed3; targets_bed4 }
Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2; ref_fasta3 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2; ref_fai3 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2; ref_dict3 }
Channel.fromPath( file(params.ref_chrom_sizes) ).set{ ref_chrom_sizes }
Channel.fromPath( file(params.trimmomatic_contaminant_fa) ).set{ trimmomatic_contaminant_fa }
Channel.fromPath( file(params.ref_fa_bwa_dir) ).set{ ref_fa_bwa_dir }
Channel.fromPath( file(params.gatk_1000G_phase1_indels_hg19_vcf) ).set{ gatk_1000G_phase1_indels_vcf }
Channel.fromPath( file(params.mills_and_1000G_gold_standard_indels_hg19_vcf) ).set{ mills_and_1000G_gold_standard_indels_vcf }
Channel.fromPath( file(params.dbsnp_ref_vcf) ).into{ dbsnp_ref_vcf; dbsnp_ref_vcf2; dbsnp_ref_vcf3 }
Channel.fromPath( file(params.cosmic_ref_vcf) ).into{ cosmic_ref_vcf; cosmic_ref_vcf2 }
Channel.fromPath( file(params.microsatellites) ).set{ microsatellites }



// read samples from analysis samplesheet
Channel.fromPath( file(params.samples_analysis_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_ID = row['Sample']
            def reads1 = row['R1'].tokenize( ',' ).collect { file(it) } // comma-sep string into list of files
            def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
            return [ sample_ID, reads1, reads2 ]
        }
        .tap { samples_R1_R2; samples_R1_R2_2 } // set of all fastq R1 R2 per sample
        .map { sample_ID, reads1, reads2 ->
            return [ reads1, reads2 ]
        }
        .flatMap().flatMap()
        .set { samples_each_fastq } // emit each fastq file individually, no sampleID


// read sample tumor-normal pairs from analysis sheet
Channel.fromPath( file(params.samples_analysis_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def tumorID = row['Tumor']
            def normalID = row['Normal']
            return [ tumorID, normalID ]
        }
        .filter { item ->
            item[0] != item[1] // remove Normal samples
        }
        .filter { item ->
            item[1] != 'NA' // unpaired samples
        }
        .into { samples_pairs; samples_pairs2 }

// view paired entries
samples_pairs2.subscribe { println "samples_pairs2: ${it}" }





//
// PIPELINE TASKS
//

// PREPROCESSING
process fastq_merge {
    // merge the R1 and R2 fastq files into a single fastq each
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/fastq-merge", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(fastq_r1: "*"), file(fastq_r2: "*") from samples_R1_R2

    output:
    set val(sample_ID), file("${sample_ID}_R1.fastq.gz"), file("${sample_ID}_R2.fastq.gz") into samples_fastq_merged

    script:
    """
    cat ${fastq_r1} > "${sample_ID}_R1.fastq.gz"
    cat ${fastq_r2} > "${sample_ID}_R2.fastq.gz"
    """
}


process trimmomatic {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/fastq-trim", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(read1), file(read2), file(trimmomatic_contaminant_fa) from samples_fastq_merged.combine(trimmomatic_contaminant_fa)

    output:
    set val(sample_ID), file("${sample_ID}_R1.trim.fastq.gz"), file("${sample_ID}_R2.trim.fastq.gz") into samples_fastq_trimmed, samples_fastq_trimmed2

    script:
    """
    trimmomatic.sh PE -threads \${NSLOTS:-\${NTHREADS:-1}} \
    "${read1}" "${read2}" \
    "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
    "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:${trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
    """
}


process fastqc_trim {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/fastqc-trim", mode: 'copy', overwrite: true

    input:
    set val(sample_ID),  file(fastq_R1_trim), file(fastq_R2_trim) from samples_fastq_trimmed2

    output:
    file(output_R1_html)
    file(output_R1_zip)
    file(output_R2_html)
    file(output_R2_zip)

    script:
    output_R1_html = "${fastq_R1_trim}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
    output_R1_zip = "${fastq_R1_trim}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
    output_R2_html = "${fastq_R2_trim}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
    output_R2_zip = "${fastq_R2_trim}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
    """
    fastqc -o . "${fastq_R1_trim}"
    fastqc -o . "${fastq_R2_trim}"
    """

}


process bwa_mem {
    // first pass alignment with BWA
    tag { "${sample_ID}" }

    input:
    set val(sample_ID), file(fastq_R1_trim), file(fastq_R2_trim), file(ref_fa_bwa_dir) from samples_fastq_trimmed.combine(ref_fa_bwa_dir)

    output:
    set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam

    script:
    """
    bwa mem -M -v 1 -t \${NSLOTS:-\${NTHREADS:-1}} -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${ref_fa_bwa_dir}/genome.fa" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
    """
}
