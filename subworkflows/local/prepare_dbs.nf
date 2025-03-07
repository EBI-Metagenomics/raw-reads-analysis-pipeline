/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
include { GET_RFAM_DB } from '../modules/local/download_dbs/rfam_models/main'
include { MOTUS_DOWNLOAD_DB } from '../modules/nf-core/motus/download_db/main'
include { GET_REFERENCE_GENOME } from '../modules/local/download_dbs/reference_genome/main'
include { GET_MAPSEQ_DB } from '../modules/local/download_dbs/mapseq_db/main'

workflow DOWNLOAD_RFAM {
    main:
        def meta = [id: params.databases.rfam.name]
        if (params.databases.rfam.local_path.size()>0) {
            local_path = "${params.databases.rfam.local_path}"
        }else{
            local_path = false
        }
        if ( local_path && file("${local_path}").exists() ) {
            rfam_db_dir = Channel.value([meta, file(local_path)])
        }else{
            cache_path = "${params.databases.cache_path}/${params.databases.rfam.name}"
            if ( file("${cache_path}").exists() ) {    
                rfam_db_dir = Channel.value([meta, file(cache_path)])
            }else{
                GET_RFAM_DB([meta, "${params.databases.rfam.remote_path}", "${params.databases.rfam.md5_path}"])
                rfam_db_dir = GET_RFAM_DB.out.db_dir.first()
            }
        }

    emit:
        rfam_db_dir
}

workflow DOWNLOAD_MOTUS_DB {
    main:
        def meta = [id: params.databases.motus.name]
        if (params.databases.motus.local_path.size()>0) {
            local_path = "${params.databases.motus.local_path}"
        }else{
            local_path = false
        }
        if ( local_path && file("${local_path}").exists() ) {
            motus_db = Channel.value([meta, file(local_path)])
        }else{
            cache_path = "${params.databases.cache_path}/${params.databases.motus.name}"
            if ( file("${cache_path}").exists() ) {    
                motus_db = Channel.value([meta, file(cache_path)])
            }else{
                // MOTUS_DOWNLOAD_DB("${params.databases.motus.remote_path}", "${params.databases.motus.md5_path}")
                MOTUS_DOWNLOAD_DB()
                motus_db = MOTUS_DOWNLOAD_DB.out.db
            }
        }
    emit:
        motus_db
}

workflow DOWNLOAD_HOST_REFERENCE_GENOME {
    main:
        def meta = [id: params.databases.host_genome.name]
        if (params.databases.host_genome.local_path.size()>0) {
            local_path = "${params.databases.host_genome.local_path}"
        }else{
            local_path = false
        }
        if ( local_path && file("${local_path}/${params.databases.host_genome.files.genome}").exists() ) {
            host_genome_db = Channel.value([meta, file(local_path)])
        }else{
            cache_path = "${params.databases.cache_path}/${params.databases.host_genome.name}"
            if ( file("${cache_path}/${params.databases.host_genome.files.genome}").exists() ) {    
                host_genome_db = Channel.value([meta, file(cache_path)])
            }else{
                GET_REFERENCE_GENOME([meta, "${params.databases.host_genome.remote_path}", "${params.databases.host_genome.md5_path}"])
                host_genome_db = GET_REFERENCE_GENOME.out.first()
            }
        }
    emit:
        host_genome_db
}

workflow DOWNLOAD_MAPSEQ_SSU {
    main:
        def meta = [id: 'silva_ssu', name: params.databases.silva_ssu.name]
        if (params.databases.silva_ssu.local_path.size()>0) {
            local_path = "${params.databases.silva_ssu.local_path}"
        }else{
            local_path = false
        }
        if ( local_path && file("${local_path}/${params.databases.silva_ssu.name}").exists() ) {
            mapseq_ssu_db_dir = Channel.value([meta, file(local_path)])
        }else{
            cache_path = "${params.databases.cache_path}/${params.databases.silva_ssu.name}"
            if ( file("${cache_path}/${params.databases.silva_ssu.name}").exists() ) {    
                mapseq_ssu_db_dir = Channel.value([meta, file(cache_path)])
            }else{
                GET_MAPSEQ_DB([
                    meta,
                    "${params.databases.silva_ssu.remote_path}", 
                    "${params.databases.silva_ssu.md5_path}"
                ])
                mapseq_ssu_db_dir = GET_MAPSEQ_DB.out.db_dir
            }
        }
    emit:
        mapseq_ssu_db_dir
}

workflow DOWNLOAD_MAPSEQ_LSU {
    main:
        def meta = [id: 'silva_lsu', name: params.databases.silva_lsu.name]
        if (params.databases.silva_lsu.local_path.size()>0) {
            local_path = "${params.databases.silva_lsu.local_path}"
        }else{
            local_path = false
        }
        if ( local_path && file("${local_path}/${params.databases.silva_lsu.name}").exists() ) {
            mapseq_lsu_db_dir = Channel.value([meta, file(local_path)])
        }else{
            cache_path = "${params.databases.cache_path}/${params.databases.silva_lsu.name}"
            if ( file("${cache_path}/${params.databases.silva_lsu.name}").exists() ) {    
                mapseq_lsu_db_dir = Channel.value([file(cache_path)])
            }else{
                GET_MAPSEQ_DB([
                    meta,
                    "${params.databases.silva_lsu.remote_path}",
                    "${params.databases.silva_lsu.md5_path}"
                ])
                mapseq_lsu_db_dir = GET_MAPSEQ_DB.out.db_dir
            }
        }
    emit:
        mapseq_lsu_db_dir
}

