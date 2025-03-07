/*
 * Host decontamination
*/
process DECONTAMINATION {

    publishDir "${params.outdir}/qc/decontamination", mode: 'copy'
    
    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'
    label 'decontamination'
    tag "$meta.id"

    input:
    tuple val(meta), val(fasta)
    tuple val(db_meta), val(host_genome_db)

    output:
    tuple val(meta), path("*_clean*.fastq.gz")

    script:
    def input_reads = "";
    def bwa_index = "${host_genome_db}/${params.databases.host_genome.files.bwa_index_prefix}"
    if (meta.mode == "single") {
        input_reads = "${fasta}";
        """
        mkdir -p output_decontamination

        echo "mapping files to host genome SE"
        bwa-mem2 mem -t ${task.cpus} \
        ${bwa_index} \
        ${input_reads} > out.sam

        echo "convert sam to bam"
        samtools view -@ ${task.cpus} -f 4 -F 256 -uS -o output_decontamination/${meta.id}_unmapped.bam out.sam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${meta.id}_unmapped.bam \
        -o output_decontamination/${meta.id}_unmapped_sorted.bam

        echo "samtools"
        samtools fastq output_decontamination/${meta.id}_unmapped_sorted.bam > output_decontamination/${meta.id}_clean.fastq

        echo "compressing output file"
        gzip -c output_decontamination/${meta.id}_clean.fastq > ${meta.id}_clean.fastq.gz
        """
    } else if ( meta.mode == "paired" ) {
        if (fasta[0].name.contains("_1")) {
            input_reads = "${fasta[0]} ${fasta[1]}"
        } else {
            input_reads = "${fasta[1]} ${fasta[0]}"
        }
        """
        mkdir output_decontamination
        echo "mapping files to host genome PE"
        bwa-mem2 mem \
        -t ${task.cpus} \
        ${bwa_index} \
        ${input_reads} > out.sam

        echo "convert sam to bam"
        samtools view -@ ${task.cpus} -f 12 -F 256 -uS -o output_decontamination/${meta.id}_both_unmapped.bam out.sam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${meta.id}_both_unmapped.bam -o output_decontamination/${meta.id}_both_unmapped_sorted.bam

        echo "samtools fastq"
        samtools fastq -1 output_decontamination/${meta.id}_clean_1.fastq \
        -2 output_decontamination/${meta.id}_clean_2.fastq \
        -0 /dev/null \
        -s /dev/null \
        -n output_decontamination/${meta.id}_both_unmapped_sorted.bam

        echo "compressing output files"
        gzip -c output_decontamination/${meta.id}_clean_1.fastq > ${meta.id}_clean_1.fastq.gz
        gzip -c output_decontamination/${meta.id}_clean_2.fastq > ${meta.id}_clean_2.fastq.gz
        """
    } else {
        error "Invalid mode: ${sample_d.mode}"
    }
}

/*
 * Decontamination counts
 * # TODO: assess Qualimap to see tha bwa mapping results
*/
process DECONTAMINATION_REPORT {

    publishDir "${params.outdir}/qc/decontamination", mode: 'copy'

    label 'decontamination_report'
    tag "$meta.id"
    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), val(cleaned_reads)

    output:
    tuple val(meta), path("decontamination_output_report.txt")

    script:
    def input_f_reads = "";
    def input_r_reads = "";
    if ( meta.mode == "paired" ) {
        if (cleaned_reads[0].name.contains('_1')) {
            input_f_reads = cleaned_reads[0]
            input_r_reads = cleaned_reads[1]
        } else if (cleaned_reads[1].name.contains('_1')) {
            input_f_reads = cleaned_reads[1]
            input_r_reads = cleaned_reads[0]
        }
        """
        zcat ${input_f_reads} | grep '@' | wc -l > decontamination_output_report.txt
        zcat ${input_r_reads} | grep '@' | wc -l >> decontamination_output_report.txt
        """
    } else if ( meta.mode == "single" ) {
        """
        zcat ${cleaned_reads} | grep '@' | wc -l > decontamination_output_report.txt
        """
    } else {
        error "Invalid mode: ${meta.mode}"
    }
}



