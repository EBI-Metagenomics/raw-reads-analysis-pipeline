include { ALIGNHOST } from '../../modules/local/alignhost/main'

workflow DECONTAMINATION {
    take:
    reads            // tuple(meta, reads)
    ref_genome       // path(reference_genome)
    ref_genome_index // path(reference_genome_index

    main:

    ch_versions = Channel.empty()

    to_align = reads.map { meta, reads -> 
        [ meta, reads, ref_genome, ref_genome_index ]
    }

    ALIGNHOST( to_align )

    ch_versions = ch_versions.mix(ALIGNHOST.out.versions.first())

    emit:
    decontaminated_reads = ALIGNHOST.out.reads
    versions = ch_versions                          // channel: [ versions.yml ]
}

