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
process fastqc_raw {
    tag { "${fastq}" }
    module "fastqc/0.11.7"
    publishDir "${params.output_dir}/fastqc-raw", mode: 'copy', overwrite: true

    input:
    file(fastq) from samples_each_fastq

    output:
    file(output_html)
    file(output_zip)

    script:
    output_html = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
    output_zip = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
    """
    echo "output_zip: ${output_zip}, output_html: ${output_html}"
    fastqc -o . "${fastq}"
    """

}

process fastq_merge {
    // merge the R1 and R2 fastq files into a single fastq each
    tag { "${sample_ID}" }
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
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
    // Illumina read trimming
    // http://www.usadellab.org/cms/?page=trimmomatic
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/fastq-trim", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8 -l mem_free=40G -l mem_token=5G'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(read1), file(read2), file(trimmomatic_contaminant_fa) from samples_fastq_merged.combine(trimmomatic_contaminant_fa)

    output:
    set val(sample_ID), file("${sample_ID}_R1.trim.fastq.gz"), file("${sample_ID}_R2.trim.fastq.gz") into samples_fastq_trimmed, samples_fastq_trimmed2

    script:
    """
    java -Xms16G -Xmx16G -jar ${params.trimmomatic_jar} PE -threads \${NSLOTS:-1} \
    "${read1}" "${read2}" \
    "${sample_ID}_R1.trim.fastq.gz" "${sample_ID}_R1.unpaired.fastq.gz" \
    "${sample_ID}_R2.trim.fastq.gz" "${sample_ID}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:${trimmomatic_contaminant_fa}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
    """
}

process fastqc_trim {
    tag { "${sample_ID}" }
    module "fastqc/0.11.7"
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
    clusterOptions '-pe threaded 4-16 -l mem_free=40G -l mem_token=4G'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'bwa/0.7.17'

    input:
    set val(sample_ID), file(fastq_R1_trim), file(fastq_R2_trim), file(ref_fa_bwa_dir) from samples_fastq_trimmed.combine(ref_fa_bwa_dir)

    output:
    set val(sample_ID), file("${sample_ID}.sam") into samples_bwa_sam

    script:
    """
    bwa mem -M -v 1 -t \${NSLOTS:-1} -R '@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA' "${ref_fa_bwa_dir}/genome.fa" "${fastq_R1_trim}" "${fastq_R2_trim}" -o "${sample_ID}.sam"
    """
}

process sambamba_view_sort {
    tag { "${sample_ID}" }
    clusterOptions '-pe threaded 1-8 -l mem_free=40G -l mem_token=4G'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_sam) from samples_bwa_sam

    output:
    set val(sample_ID), file("${sample_ID}.bam") into samples_bam, samples_bam2

    script:
    """
    "${params.sambamba_bin}" view --sam-input --nthreads=\${NSLOTS:-1} --filter='mapping_quality>=10' --format=bam --compression-level=0 "${sample_sam}" | \
    "${params.sambamba_bin}" sort --nthreads=\${NSLOTS:-1} --memory-limit=16GB --out="${sample_ID}.bam" /dev/stdin
    """
}

process sambamba_flagstat {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/sambamba-flagstat", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_bam

    output:
    file "${sample_ID}.flagstat.txt"

    script:
    """
    "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.flagstat.txt"
    """
}

process sambamba_dedup {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/bam-bwa-dd", mode: 'copy', overwrite: true
    clusterOptions '-pe threaded 1-8 -l mem_free=40G -l mem_token=4G'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam) from samples_bam2

    output:
    set val(sample_ID), file("${sample_ID}.dd.bam") into samples_dd_bam, samples_dd_bam2, samples_dd_bam3, samples_dd_bam4, samples_dd_bam5, samples_dd_bam6, samples_dd_bam7
    file("${sample_ID}.dd.bam.bai")

    script:
    """
    "${params.sambamba_bin}" markdup --remove-duplicates --nthreads \${NSLOTS:-1} --hash-table-size 525000 --overflow-list-size 525000 "${sample_bam}" "${sample_ID}.dd.bam"
    samtools view "${sample_ID}.dd.bam"
    """
}

process sambamba_dedup_flagstat {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/sambamba-dd-flagstat", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(sample_ID), file(sample_bam) from samples_dd_bam2

    output:
    file "${sample_ID}.dd.flagstat.txt"

    script:
    """
    "${params.sambamba_bin}" flagstat "${sample_bam}" > "${sample_ID}.dd.flagstat.txt"
    """

}






// setup downstream Channels
samples_dd_bam.combine(ref_fasta)
            .combine(ref_fai)
            .combine(ref_dict)
            .tap { samples_dd_bam_ref }
            .combine(targets_bed)
            .tap { samples_dd_bam_ref2;
                    samples_dd_bam_ref3;
                    samples_dd_bam_ref4
                }
            .combine(gatk_1000G_phase1_indels_vcf)
            .combine(mills_and_1000G_gold_standard_indels_vcf)
            .combine(dbsnp_ref_vcf)
            .set { samples_dd_bam_ref_gatk }






// GATK RECALIBRATION AND VARIANT CALLING

process qc_target_reads_gatk_genome {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-16 -l mem_free=40G -l mem_token=5G'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict) from samples_dd_bam_ref

    output:
    file "${sample_ID}.genome.sample_statistics"
    file "${sample_ID}.genome.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.genome"
    """
}


process qc_target_reads_gatk_pad500 {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-16 -l mem_free=40G -l mem_token=5G'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref2

    output:
    file "${sample_ID}.pad500.sample_statistics"
    file "${sample_ID}.pad500.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 500 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.pad500"
    """
}

process qc_target_reads_gatk_pad100 {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-16 -l mem_free=40G -l mem_token=5G'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref3

    output:
    file "${sample_ID}.pad100.sample_statistics"
    file "${sample_ID}.pad100.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 100 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.pad100"
    """
}

process qc_target_reads_gatk_bed {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc-target-reads", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-16 -l mem_free=40G -l mem_token=5G'

    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_bam_ref4

    output:
    file "${sample_ID}.bed.sample_statistics"
    file "${sample_ID}.bed.sample_summary"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    -rf BadCigar \
    -nt \${NSLOTS:-1} \
    --logging_level ERROR \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 100 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}.bed"
    """
}


// MAIN REALIGNMENT AND RECALIBRATION STEP
process bam_ra_rc_gatk {
    // re-alignment and recalibration with GATK
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/bam_dd_ra_rc_gatk", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 4-16 -l mem_free=40G -l mem_token=4G'
    module 'samtools/1.3'


    input:
    set val(sample_ID), file(sample_bam), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(gatk_1000G_phase1_indels_vcf), file(mills_and_1000G_gold_standard_indels_vcf), file(dbsnp_ref_vcf) from samples_dd_bam_ref_gatk

    output:
    set val(sample_ID), file("${sample_ID}.dd.ra.rc.bam"), file("${sample_ID}.dd.ra.rc.bam.bai") into samples_dd_ra_rc_bam, samples_dd_ra_rc_bam2, samples_dd_ra_rc_bam3
    file "${sample_ID}.intervals"
    file "${sample_ID}.table1.txt"
    file "${sample_ID}.table2.txt"
    file "${sample_ID}.csv"
    file "${sample_ID}.pdf"

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T RealignerTargetCreator \
    -dt NONE \
    --logging_level ERROR \
    -nt \${NSLOTS:-1} \
    --reference_sequence "${ref_fasta}" \
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.intervals"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T IndelRealigner \
    -dt NONE \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    --maxReadsForRealignment 50000 \
    -known "${gatk_1000G_phase1_indels_vcf}" \
    -known "${mills_and_1000G_gold_standard_indels_vcf}" \
    -targetIntervals "${sample_ID}.intervals" \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.dd.ra.bam"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.table1.txt"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T BaseRecalibrator \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -knownSites "${gatk_1000G_phase1_indels_vcf}" \
    -knownSites "${mills_and_1000G_gold_standard_indels_vcf}" \
    -knownSites "${dbsnp_ref_vcf}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_ID}.dd.ra.bam" \
    -BQSR "${sample_ID}.table1.txt" \
    --out "${sample_ID}.table2.txt"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T AnalyzeCovariates \
    --logging_level ERROR \
    --reference_sequence "${ref_fasta}" \
    -before "${sample_ID}.table1.txt" \
    -after "${sample_ID}.table2.txt" \
    -csv "${sample_ID}.csv" \
    -plots "${sample_ID}.pdf"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T PrintReads \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    -BQSR "${sample_ID}.table1.txt" \
    --input_file "${sample_ID}.dd.ra.bam" \
    --out "${sample_ID}.dd.ra.rc.bam"

    samtools index "${sample_ID}.dd.ra.rc.bam"
    """
}







// setup downstream Channels for per-sample analyses
samples_dd_ra_rc_bam.combine(ref_fasta2)
                    .combine(ref_fai2)
                    .combine(ref_dict2)
                    .combine(targets_bed2)
                    .tap { samples_dd_ra_rc_bam_ref;
                            samples_dd_ra_rc_bam_ref2;
                            samples_dd_ra_rc_bam_ref3;
                            samples_dd_ra_rc_bam_ref4;
                            samples_dd_ra_rc_bam_ref5;
                            samples_dd_ra_rc_bam_ref6;
                            samples_dd_ra_rc_bam_ref7;
                            samples_dd_ra_rc_bam_ref8 }


samples_dd_ra_rc_bam_ref2.combine( dbsnp_ref_vcf3 )
                        .tap { samples_dd_ra_rc_bam_ref_dbsnp;
                                samples_dd_ra_rc_bam_ref_dbsnp2 }









process qc_coverage_gatk {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/qc_coverage_gatk", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 1-16 -l mem_free=40G -l mem_token=5G'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref

    output:
    file "${sample_ID}.sample_summary"
    file "${sample_ID}.sample_statistics"
    file "${sample_ID}.sample_interval_summary"
    file "${sample_ID}.sample_interval_statistics"
    file "${sample_ID}.sample_cumulative_coverage_proportions"
    file "${sample_ID}.sample_cumulative_coverage_counts"
    file("${sample_ID}.summary.csv") into qc_coverage_gatk_summary

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T DepthOfCoverage \
    -dt NONE \
    --logging_level ERROR \
    -rf BadCigar \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --omitDepthOutputAtEachBase \
    -ct 10 -ct 50 -ct 100 -ct 500 \
    --minBaseQuality 20 \
    --minMappingQuality 20 \
    --nBins 999 \
    --start 1 --stop 1000 \
    --input_file "${sample_bam}" \
    --outputFormat csv \
    --out "${sample_ID}"

    head -2 "${sample_ID}.sample_summary" > "${sample_ID}.summary.csv"
    """
}
qc_coverage_gatk_summary.collectFile(name: "${params.qc_coverage_gatk_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)

process pad_bed {
    publishDir "${params.output_dir}/targets", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    module 'bedtools/2.26.0'

    input:
    set file(targets_bed_file), file(ref_chrom_sizes) from targets_bed3.combine(ref_chrom_sizes)

    output:
    file("targets.pad10.bed") into targets_pad_bed

    script:
    """
    cat "${targets_bed_file}" | LC_ALL=C sort -k1,1 -k2,2n | bedtools slop -g "${ref_chrom_sizes}" -b 10 | bedtools merge -d 5 > targets.pad10.bed
    """
}

process lofreq {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/vcf_lofreq", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 4-16 -l mem_free=40G -l mem_token=4G'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(dbsnp_ref_vcf) from samples_dd_ra_rc_bam_ref_dbsnp

    output:
    file("${sample_ID}.vcf")
    file("${sample_ID}.norm.vcf")
    file("${sample_ID}.norm.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt") into lofreq_annotations
    file("${sample_ID}.eval.grp")
    val(sample_ID) into sample_lofreq_done

    script:
    """
    "${params.lofreq_bin}" call-parallel \
    --call-indels \
    --pp-threads \${NSLOTS:-1} \
    --ref "${ref_fasta}" \
    --bed "${targets_bed_file}" \
    --out "${sample_ID}.vcf" \
    "${sample_bam}"

    bgzip -c "${sample_ID}.vcf" > "${sample_ID}.vcf.bgz"

    bcftools index "${sample_ID}.vcf.bgz"

    bcftools norm \
    --multiallelics \
    -both \
    --output-type v \
    "${sample_ID}.vcf.bgz" | \
    bcftools norm \
    --fasta-ref "${ref_fasta}" \
    --output-type v - | \
    bcftools view \
    --exclude 'DP<5' \
    --output-type v >  "${sample_ID}.norm.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.norm.vcf" "${sample_ID}.norm"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.norm.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.norm.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T VariantEval \
    -R "${ref_fasta}" \
    -o "${sample_ID}.eval.grp" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --eval "${sample_ID}.norm.vcf"

    """
}
lofreq_annotations.collectFile(name: "${params.annotations_lofreq_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)




process gatk_hc {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/vcf_hc", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-pe threaded 4-16 -l mem_free=40G -l mem_token=4G'
    module 'samtools/1.3'

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file), file(dbsnp_ref_vcf) from samples_dd_ra_rc_bam_ref_dbsnp2

    output:
    file("${sample_ID}.vcf")
    set val(sample_ID), file("${sample_ID}.norm.vcf") into sample_vcf_hc
    file("${sample_ID}.norm.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt") into gatk_hc_annotations
    set val(sample_ID), file("${sample_ID}.norm.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt") into gatk_hc_annotations2
    file("${sample_ID}.eval.grp")
    val(sample_ID) into sample_gatk_hc_done

    script:
    """
    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T HaplotypeCaller \
    -dt NONE \
    --logging_level ERROR \
    -nct \${NSLOTS:-1} \
    --max_alternate_alleles 3 \
    --standard_min_confidence_threshold_for_calling 50 \
    --reference_sequence "${ref_fasta}" \
    --intervals "${targets_bed_file}" \
    --interval_padding 10 \
    --input_file "${sample_bam}" \
    --out "${sample_ID}.vcf"

    cat "${sample_ID}.vcf" | \
    bcftools norm \
    --multiallelics \
    -both \
    --output-type v - | \
    bcftools norm \
    --fasta-ref "${ref_fasta}" \
    --output-type v - | \
    bcftools view \
    --exclude 'DP<5' \
    --output-type v > "${sample_ID}.norm.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.norm.vcf" "${sample_ID}.norm"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.norm.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.norm.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T VariantEval \
    -R "${ref_fasta}" \
    -o "${sample_ID}.eval.grp" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --eval "${sample_ID}.norm.vcf"
    """
}
gatk_hc_annotations.collectFile(name: "${params.annotations_hc_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)






//
// DOWNSTREAM TASKS
//
// DELLY2 SNV STEPS
process delly2_deletions {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/snv-deletions-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref4

    output:
    file "${sample_ID}.deletions.vcf"
    file "${sample_ID}.deletions.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" into delly2_deletions_annotations

    script:
    """
    "${params.delly2_bin}" call -t DEL -g "${ref_fasta}" -o "${sample_ID}.deletions.bcf" "${sample_bam}"
    "${params.delly2_bcftools_bin}" view "${sample_ID}.deletions.bcf" > "${sample_ID}.deletions.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.deletions.vcf" "${sample_ID}.deletions"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.deletions.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.deletions.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}
delly2_deletions_annotations.collectFile(name: "${params.annotations_deletions_Delly2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)


process delly2_duplications {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/snv-duplications-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref5

    output:
    file "${sample_ID}.duplications.vcf"
    file "${sample_ID}.duplications.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" into delly2_duplications_annotations

    script:
    """
    "${params.delly2_bin}" call -t DUP -g "${ref_fasta}" -o "${sample_ID}.duplications.bcf" "${sample_bam}"
    "${params.delly2_bcftools_bin}" view "${sample_ID}.duplications.bcf" > "${sample_ID}.duplications.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.duplications.vcf" "${sample_ID}.duplications"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.duplications.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.duplications.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}
delly2_duplications_annotations.collectFile(name: "${params.annotations_duplications_Delly2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)



process delly2_inversions {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/snv-inversions-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref6

    output:
    file "${sample_ID}.inversions.bcf"
    file "${sample_ID}.inversions.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" into delly2_inversions_annotations

    script:
    """
    "${params.delly2_bin}" call -t INV -g "${ref_fasta}" -o "${sample_ID}.inversions.bcf" "${sample_bam}"
    "${params.delly2_bcftools_bin}" view "${sample_ID}.inversions.bcf" > "${sample_ID}.inversions.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.inversions.vcf" "${sample_ID}.inversions"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.inversions.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.inversions.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}
delly2_inversions_annotations.collectFile(name: "${params.annotations_inversion_Delly2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)



process delly2_translocations {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/snv-translocations-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref7

    output:
    file "${sample_ID}.translocations.vcf"
    file "${sample_ID}.translocations.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" into delly2_translocations_annotations

    script:
    """
    ${params.delly2_bin} call -t BND -g ${ref_fasta} -o "${sample_ID}.translocations.bcf" "${sample_bam}"
    ${params.delly2_bcftools_bin} view "${sample_ID}.translocations.bcf" > "${sample_ID}.translocations.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.translocations.vcf" "${sample_ID}.translocations"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.translocations.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.translocations.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}
delly2_translocations_annotations.collectFile(name: "${params.annotations_translocations_Delly2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)


process delly2_insertions {
    tag { "${sample_ID}" }
    publishDir "${params.output_dir}/snv-insertions-Delly2", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_bam), file(sample_bai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed_file) from samples_dd_ra_rc_bam_ref8

    output:
    file "${sample_ID}.insertions.vcf"
    file "${sample_ID}.insertions.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" into delly2_insertions_annotations

    script:
    """
    ${params.delly2_bin} call -t INS -g ${ref_fasta} -o "${sample_ID}.insertions.bcf" "${sample_bam}"
    ${params.delly2_bcftools_bin} view "${sample_ID}.insertions.bcf" > "${sample_ID}.insertions.vcf"

    # annotate the vcf
    annotate_vcf.sh "${sample_ID}.insertions.vcf" "${sample_ID}.insertions"

    # add a column with the sample ID
    paste_col.py -i "${sample_ID}.insertions.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${sample_ID}.insertions.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${sample_ID}" -d "\t"
    """
}
delly2_insertions_annotations.collectFile(name: "${params.annotations_insertions_Delly2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)



// Genomic Signatures
process deconstructSigs_signatures {
    tag { "${sample_ID}" }
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    publishDir "${params.output_dir}/signatures_hc", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_vcf) from sample_vcf_hc

    output:
    file "${sample_ID}_signatures.Rds"
    file "${sample_ID}_signatures.pdf" into signatures_plots
    file "${sample_ID}_signatures_pie.pdf" into signatures_pie_plots

    script:
    """
    deconstructSigs_make_signatures.R "${sample_ID}" "${sample_vcf}"
    """
}


process merge_signatures_plots {
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'
    executor "local"
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    file '*' from signatures_plots.toList()

    output:
    file "signatures.pdf"

    script:
    """
    if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
        gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=genomic_signatures.pdf *
    else
        exit 11
    fi
    """
}


process merge_signatures_pie_plots {
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'
    executor "local"
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    file '*' from signatures_pie_plots.toList()

    output:
    file "signatures_pie.pdf"

    script:
    """
    if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=genomic_signatures_pie.pdf *
    else
        exit 11
    fi
    """
}


// // REQUIRES ANNOTATIONS FOR DBSNP FILTERING
process vaf_distribution_plot {
    tag { "${sample_ID}" }
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'
    publishDir "${params.output_dir}/vaf-distribution-hc", mode: 'copy', overwrite: true

    input:
    set val(sample_ID), file(sample_vcf_annot) from gatk_hc_annotations2

    output:
    file "${sample_ID}_vaf_dist.pdf" into vaf_distribution_plots

    script:
    """
    VAF-distribution-plot.R "${sample_ID}" "${sample_vcf_annot}"
    """

}

process merge_VAF_plots {
    executor "local"
    validExitStatus 0,11 // allow '11' failure triggered by few/no variants
    errorStrategy 'ignore'
    publishDir "${params.output_dir}/", mode: 'copy', overwrite: true

    input:
    file '*' from vaf_distribution_plots.toList()

    output:
    file "vaf_distributions.pdf"

    script:
    """
    if [ "\$(ls -1 * | wc -l)" -gt 0 ]; then
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=vaf_distributions.pdf *
    else
        exit 11
    fi
    """
}




// SETUP CHANNELS FOR PAIRED TUMOR-NORMAL STEPS
// samples_dd_ra_rc_bam2 // set val(sample_ID), file("${sample_ID}.dd.ra.rc.bam"), file("${sample_ID}.dd.ra.rc.bam.bai")
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
                    .tap { samples_dd_ra_rc_bam_pairs_ref_msi }




// get the unique chromosomes in the targets bed file
//  for per-chrom paired variant calling
Channel.fromPath( params.targets_bed )
            .splitCsv(sep: '\t')
            .map{row ->
                row[0]
            }
            .unique()
            .set{ chroms }

// add the reference .vcf
samples_dd_ra_rc_bam_pairs_ref.combine(dbsnp_ref_vcf2)
                            .combine(cosmic_ref_vcf2)
                            .tap { samples_dd_ra_rc_bam_pairs_ref_gatk } // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), dbsnp, cosmic ]
                            // add the chroms
                            .combine( chroms ) // [ comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), dbsnp, cosmic, chrom ]
                            .set { samples_dd_ra_rc_bam_pairs_ref_gatk_chrom }



samples_dd_ra_rc_bam_pairs2.subscribe { println "samples_dd_ra_rc_bam_pairs2: ${it}" }

process tumor_normal_compare {
    tag { "${comparisonID}" }
    echo true
    executor "local"
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed) from samples_dd_ra_rc_bam_pairs_ref3

    script:
    """
    echo "[tumor_normal_compare] comparisonID: ${comparisonID}, tumorID: ${tumorID}, tumorBam: ${tumorBam}, tumorBai: ${tumorBai}, normalID: ${normalID}, normalBam: ${normalBam}, normalBai: ${normalBai}, ref_fasta: ${ref_fasta}, ref_fai: ${ref_fai}, ref_dict: ${ref_dict}, targets_bed: ${targets_bed}, "
    """
}


// REQUIRES PAIRED SAMPLES BAM FILES
process msisensor {
    tag { "${comparisonID}" }
    module 'samtools/1.3'
    clusterOptions '-pe threaded 1-8 -l mem_free=40G -l mem_token=5G'
    publishDir "${params.output_dir}/microsatellites", mode: 'copy', overwrite: true


    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(microsatellites) from samples_dd_ra_rc_bam_pairs_ref_msi

    output:
    file "${comparisonID}.msisensor"
    file "${comparisonID}.msisensor_dis"
    file "${comparisonID}.msisensor_germline"
    file "${comparisonID}.msisensor_somatic"

    script:
    """
    "${params.msisensor_bin}" msi -d "${microsatellites}" -n "${normalBam}" -t "${tumorBam}" -e "${targets_bed}" -o "${comparisonID}.msisensor" -l 1 -q 1 -b \${NSLOTS:-1}
    """
}


process mutect2 {
    tag { "${comparisonID}:${chrom}" }
    publishDir "${params.output_dir}/vcf_mutect2", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    clusterOptions '-l mem_free=150G -hard'
    module 'samtools/1.3'
    module 'java/1.8'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(dbsnp_ref_vcf), file(cosmic_ref_vcf), val(chrom) from samples_dd_ra_rc_bam_pairs_ref_gatk_chrom

    output:
    file("${comparisonID}.${chrom}.vcf")
    file("${comparisonID}.${chrom}.sample.chrom.${params.ANNOVAR_BUILD_VERSION}_multianno.txt") into mutect2_annotations
    val(comparisonID) into mutect2_sampleIDs

    script:
    """
    subset_bed.py "${chrom}" "${targets_bed}" > "${comparisonID}.${chrom}.bed"

    java -Xms16G -Xmx16G -jar "${params.gatk_bin}" -T MuTect2 \
    -dt NONE \
    --logging_level WARN \
    --standard_min_confidence_threshold_for_calling 30 \
    --max_alt_alleles_in_normal_count 10 \
    --max_alt_allele_in_normal_fraction 0.05 \
    --max_alt_alleles_in_normal_qscore_sum 40 \
    --reference_sequence "${ref_fasta}" \
    --dbsnp "${dbsnp_ref_vcf}" \
    --cosmic "${cosmic_ref_vcf}" \
    --intervals "${comparisonID}.${chrom}.bed" \
    --interval_padding 10 \
    --input_file:tumor "${tumorBam}" \
    --input_file:normal "${normalBam}" \
    --out "${comparisonID}.${chrom}.vcf"

    # annotate the vcf
    annotate_vcf.sh "${comparisonID}.${chrom}.vcf" "${comparisonID}.${chrom}"

    # add a column with the sample ID
    paste_col.py -i "${comparisonID}.${chrom}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${comparisonID}.${chrom}.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "Sample" -v "${comparisonID}" -d "\t"

    # add the col for this chrom
    paste_col.py -i "${comparisonID}.${chrom}.sample.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" -o "${comparisonID}.${chrom}.sample.chrom.${params.ANNOVAR_BUILD_VERSION}_multianno.txt" --header "SampleChrom" -v "${chrom}" -d "\t"
    """
}
mutect2_annotations.collectFile(name: "${params.annotations_mutect2_file_basename}", storeDir: "${params.output_dir}", keepHeader: true)


process multiqc {
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    beforeScript "${params.beforeScript_str}"
    afterScript "${params.afterScript_str}"
    executor "local"
    module 'python/2.7.3'

    input:
    val(comparisonID) from mutect2_sampleIDs.mix(sample_gatk_hc_done)
                                            .mix(sample_lofreq_done)
                                            .collect() // force it to wait for all steps to finish
    file(output_dir) from Channel.fromPath("${params.output_dir}")

    output:
    file "multiqc_report.html" into email_files
    file "multiqc_data"

    script:
    """
    export PS=\${PS:-''} # needed for virtualenv bug
    export PS1=\${PS1:-''}
    unset PYTHONPATH
    source activate
    multiqc "${output_dir}"
    """
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
        Total CPU-Hours   : ${workflow.stats.getComputeTimeString() ?: '-'}
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
        Container engine  : ${workflow.containerEngine?:'-'}
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
            // files from process channels
            attach email_files.toList().getVal()
            // files from collectFile
            // attach ["${params.output_dir}/${params.annotations_mutect2_file_basename}",
            //         "${params.output_dir}/${params.qc_coverage_gatk_file_basename}",
            //         "${params.output_dir}/${params.annotations_hc_file_basename}",
            //         "${params.output_dir}/${params.annotations_lofreq_file_basename}"]

            subject "[${params.workflow_label}] Pipeline Completion: ${status}"

            body
            """
            ${msg}
            """
            .stripIndent()
        }
    }
}
