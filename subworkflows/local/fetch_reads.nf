// include { FETCHTOOL_RAWREADS } from '../modules/fetchtool'
include { FETCH_URL } from '../modules/fetchurl'

workflow FETCH_READS {
    take:
        fetch_ch

    main:
        // check if file is remote
        remote_ch = fetch_ch
            .filter{ meta,fp,md5,t -> fp =~ /^[a-zA-Z]{2,}:\/\// }
        local_ch = fetch_ch
            .filter{ meta,fp,mg5,t -> !(fp =~ /^[a-zA-Z]{2,}:\/\//) }
            .map{ meta,fp,md5,t -> [meta,fp,t] }
        cache_path_ch = remote_ch.map{ 
            meta,fp,md5,t -> 
            [
                meta, fp, md5, t,
                file("${params.reads_cache_path}/${meta.id}/${fp.tokenize('/').last()}")
            ] 
        }
        download_ch = cache_path_ch
            .filter{ meta,fp,md5,t,local_fp -> !local_fp.exists() }
            .map{ meta,fp,md5,t,local_fp -> tuple(meta,fp,md5,t) }
        cache_ch = cache_path_ch
            .filter{ meta,fp,md5,t,local_fp -> local_fp.exists() }
            .map{ meta,fp,md5,t,local_fp -> tuple(meta,local_fp,t) } 
        downloaded_ch = FETCH_URL(download_ch)
    emit:
        // join downloaded
        local_ch.mix(cache_ch).mix(downloaded_ch)
        
}

        
