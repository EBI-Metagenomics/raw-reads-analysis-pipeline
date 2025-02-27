/*
 * Quality control, host-decontamination and overlap reads
 */

include { FASTP } from '../modules/fastp'
include { SEQPREP } from '../modules/seqprep'
include { SEQPREP_REPORT } from '../modules/seqprep'
include { QC_REPORT } from '../modules/create_qc_report'
include { DECONTAMINATION } from '../modules/decontamination'
include { DECONTAMINATION_REPORT } from '../modules/decontamination'
include { SEQTK as FASTQ_TO_FASTA} from '../modules/seqtk'
include { QC_STATS } from '../modules/qc_summary'

workflow QC {
    take:
        samples
        host_genome_db
    
    main:
        
        FASTP(samples, Channel.value(""))  // empty channel is optional merged reads
        
        fastp_out = FASTP.out.reads.map{ meta,fastq,json,html -> 
            [meta,[fastq:fastq,json:json,html:html]] 
        }
        
        DECONTAMINATION(fastp_out.join(samples), host_genome_db)
        
        decontam_out = DECONTAMINATION.out.join(samples)
        
        DECONTAMINATION_REPORT(decontam_out)

        // branch here
        seqprep_in = decontam_out.branch{ 
            meta, decontam_d, sample_d -> 
            paired: sample_d.mode=='paired'
            single: sample_d.mode=='single'
        }
        
        // merge paired reads
        SEQPREP(seqprep_in.paired)
        seqprep_out = SEQPREP.out.map{ 
            meta, merged, forward_unmapped, reverse_unmapped -> 
            [meta, [merged: merged, 
                    forward_unmapped: forward_unmapped, 
                    reverse_unmapped: reverse_unmapped]] 
        }
        SEQPREP_REPORT(seqprep_out)


        // mix single and paired samples
        paired_merged_counts = SEQPREP_REPORT.out
        single_merged_counts = seqprep_in.single
            .map{ meta, decontam_d, sample_d -> meta }
            .merge(Channel.fromPath("NO_FILE"))
        merged_counts = paired_merged_counts.mix(single_merged_counts)
        
        paired_merged_reads = seqprep_out.map{ meta,d -> [meta, d.merged] }
        single_merged_reads = seqprep_in.single.map{
            meta, decontam_d, sample_d -> [meta, decontam_d] 
        }
        merged_reads = paired_merged_reads.mix(single_merged_reads)
        
        // paired_overlapped_reads_ch = SEQPREP.out.sample_name.merge(SEQPREP.out.overlapped_reads).map{ tuple(it[0],it[1..-1]) }
        // paired_overlapped_counts_ch = SEQPREP_REPORT.out.sample_name.merge(SEQPREP_REPORT.out.overlapped_report)
        // single_overlapped_reads_ch = seqrep_ch.single.map{ tuple(it[0],it[1]) }
        // single_overlapped_counts_ch = seqrep_ch.single.map{ it[0] }.merge(Channel.fromPath("NO_FILE"))
        // overlapped_counts = paired_overlapped_counts_ch.mix(single_overlapped_counts_ch)
        // overlapped_reads = paired_overlapped_reads_ch.mix(single_overlapped_reads_ch)

        // gather reports
        fastp_report = fastp_out.map{ meta,d -> [meta,d.json] }
        decontam_report = DECONTAMINATION_REPORT.out
        report_in = merged_counts
            .join(fastp_report)
            .join(decontam_report)
            .join(samples)
            .map{ meta, a, b, c, d -> [meta, [merged: a, fastp: b, 
                                              decontam: c, samples: d]] }
        
        QC_REPORT(report_in)
        
        // generate output sequence files
        FASTQ_TO_FASTA(merged_reads)
        QC_STATS(FASTQ_TO_FASTA.out)

    emit:
        merged_reads = merged_reads
        sequence = FASTQ_TO_FASTA.out 
        qc_report = QC_REPORT.out
        qc_stats = QC_STATS.out
        fastp_json = fastp_report
}
