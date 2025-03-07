// include { FETCHTOOL_RAWREADS } from '../modules/fetchtool'
include { FETCHURL } from '../../../modules/local/fetchurl/main'

workflow FETCHREADS {
    take:
        fetch_ch  // [meta,fastq_fp]
        cache_path  // value channel
    main:
        // check if file is remote
        remote_ch = fetch_ch
            .filter{ meta,fp -> fp =~ /^[a-zA-Z]{2,}:\/\// }
        local_ch = fetch_ch
            .filter{ meta,fp -> !(fp =~ /^[a-zA-Z]{2,}:\/\//) }
            .map{ meta,fp -> [meta,fp] }
        cache_path_ch = remote_ch.map{ 
            meta,fp -> 
            [
                meta, fp,
                file("${cache_path}/${meta.id}/${fp.tokenize('/').last()}")
            ] 
        }
        download_ch = cache_path_ch
            .filter{ meta,fp,local_fp -> !local_fp.exists() }
            .map{ meta,fp,local_fp -> tuple(meta,fp) }
        cache_ch = cache_path_ch
            .filter{ meta,fp,local_fp -> local_fp.exists() }
            .map{ meta,fp,local_fp -> tuple(meta,local_fp) } 
        downloaded_ch = FETCHURL(download_ch)
    emit:
        // join downloaded
        local_ch.mix(cache_ch).mix(downloaded_ch)
        
}

        
