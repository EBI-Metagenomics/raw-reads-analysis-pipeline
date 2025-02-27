/*
 * Classify with mapseq, convert mapseq2biom, generate krona plots
 */

include { MAPSEQ } from '../modules/mapseq'
include { MAPSEQ2BIOM } from '../modules/mapseq2biom'
include { KRONA } from '../modules/krona'

workflow MAPSEQ_OTU_KRONA {
    take:
        mapseq_in  // tuple sequence, database directory

    main:
        MAPSEQ(mapseq_in)
        MAPSEQ2BIOM(MAPSEQ.out.join(mapseq_in))
        KRONA(MAPSEQ2BIOM.out)
        
        output_ch = MAPSEQ.out.join(MAPSEQ2BIOM.out)
        
    emit:
        output_ch
}

