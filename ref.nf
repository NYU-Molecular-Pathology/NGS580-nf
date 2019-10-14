Channel.from("foo").set { input_channel }

process make_ref {
    echo true
    storeDir "${params.ref_dir}"

    input:
    val(foo) from input_channel

    output:
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict")
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")

    file("BWA/hg19")
    file("Illumina/hg19/chrom.sizes")
    file("msisensor/hg19/microsatellites.list")
    file("contaminants/trimmomatic.fa")
    file("gatk-bundle/1000G_phase1.indels.hg19.vcf")
    file("gatk-bundle/1000G_phase1.indels.hg19.vcf.idx")
    file("gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf")
    file("gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf.idx")
    file("gatk-bundle/dbsnp_138.hg19.vcf")
    file("gatk-bundle/dbsnp_138.hg19.vcf.idx")
    file("hg19/CosmicCodingMuts_v73.hg19.vcf")
    file("iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")

    script:
    """
    wget https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref.tar.gz
	tar -xvzf ref.tar.gz
	rm -f ref.tar.gz
    mv ref/* .
    """
}
