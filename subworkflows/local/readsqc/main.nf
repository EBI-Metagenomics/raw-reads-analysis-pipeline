include { SEQFU_CHECK            } from "${projectDir}/modules/ebi-metagenomics/seqfu/check/main"
include { ASSESSMCPPROPORTIONS   } from "${projectDir}/modules/ebi-metagenomics/assessmcpproportions/main"
include { FASTQSUFFIXHEADERCHECK } from "${projectDir}/modules/ebi-metagenomics/fastqsuffixheadercheck/main"

workflow  READSQC {

    take:
    ch_reads    // channel: [ val(meta), [ fastq ] ]
    filter_amplicon // channel: val(boolean)

    main:
    ch_versions = Channel.empty()

    SEQFU_CHECK(ch_reads)
    ch_versions = ch_versions.mix(SEQFU_CHECK.out.versions.first())

    passed_seqfu_reads = SEQFU_CHECK.out.tsv
        .splitCsv(sep: "\t", elem: 1)
        .filter { meta, seqfu_res ->
            seqfu_res[0] == "OK"
        }
        .map { map, seqfu_res -> map }
        .join(ch_reads)

    FASTQSUFFIXHEADERCHECK(passed_seqfu_reads)
    ch_versions = ch_versions.mix(FASTQSUFFIXHEADERCHECK.out.versions.first())

    passed_suffixheader_reads = FASTQSUFFIXHEADERCHECK.out.json
        .filter { meta, sufhd_res ->
            sufhd_res.countLines() == 0
        }
        .map { meta, _ -> [ meta ] }
        .join(ch_reads)
   
    
    if ( filter_amplicon ) {
        assess_mcp_proportions_input = passed_suffixheader_reads
                        .map { meta, fastq ->
                            if ( meta.single_end ) {
                                [ meta, "auto", "auto", fastq ]
                            } else {
                                [ meta, "auto", "auto", fastq[0] ]
                            }
                        }

        ASSESSMCPPROPORTIONS(assess_mcp_proportions_input, true)
        ch_versions = ch_versions.mix(ASSESSMCPPROPORTIONS.out.versions.first())

        fastp_out = ASSESSMCPPROPORTIONS.out.library_check_out
                        .filter { meta, strategy ->
                            strategy == "AMPLICON"
                        }
                        .map { meta, _ -> [ meta ] }
                        .join(ch_reads)

        amplicon_check = ASSESSMCPPROPORTIONS.out.library_check_out
    } else {
        fastp_out = passed_suffixheader_reads
        amplicon_check = Channel.empty()
    }
 
    emit:
    passed_reads = fastp_out
    amplicon_check = amplicon_check
    versions = ch_versions               // channel: [ versions.yml ]
}
