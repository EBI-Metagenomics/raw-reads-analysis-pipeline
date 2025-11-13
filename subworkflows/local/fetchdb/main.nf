include { FETCHUNZIP } from '../../../modules/local/fetchunzip/main'

workflow FETCHDB {
    take:
    fetch_ch // channel: val(meta)
    cache_path // String

    main:
    // DBs with local paths
    local_ch = fetch_ch
        .filter { meta -> meta.local_path }
        .map { meta -> [meta, file(meta.local_path, checkIfExists: true)] }

    // if not local path and there is a remote path then the files may be in the cache
    remote_ch = fetch_ch
        .filter { meta -> ((!meta.local_path) && meta.remote_path) }
        .map { meta -> [meta, file("${cache_path}/${meta.id}")] }

    // if the remote DBs are not in the cache, or if force_download_dbs param is true then download the files 
    download_ch = remote_ch
        .filter { _meta, cache_fp -> (cache_path.isEmpty() || (!cache_fp.exists()) || params.force_download_dbs==true) }
        .map {  // Add checksum if it exists
            meta, _cache_fp -> 
            // checksum might be at apath to download from remote
            def checksum_path = []  // by default there is no checksum path
            if (meta.containsKey('remote_checksum_path')) {
                // if checksum is simply `true` then assume that the path is the remote path + '.md5' 
                if (meta.remote_checksum_path==true) {
                    checksum_path = file(meta.remote_path + '.md5', checkIfExists: true)
                } else {
                    // otherwise, if the checksum path exists then set it
                    if (meta.remote_checksum_path!=false) {
                        checksum_path = file(meta.remote_checksum_path, checkIfExists: true)
                    }
                }
            }
            // checksum can also be directly specified in parameters 
            def checksum = false
            if (meta.containsKey('remote_checksum')) {
                checksum = meta.remote_checksum
            }
            [meta, meta.id, file(meta.remote_path, checkIfExists: true), checksum_path, checksum] 
        }

    // collect channel of remote DBs that are in the cache
    cache_ch = remote_ch
       .filter { _meta, cache_fp -> ((!cache_path.isEmpty()) && cache_fp.exists() && (params.force_download_dbs==false)) }
       .map { meta, cache_fp -> [meta, cache_fp] }

    FETCHUNZIP(download_ch)
    downloaded_ch = FETCHUNZIP.out

    emit:
    dbs = local_ch.mix(cache_ch).mix(downloaded_ch)
}
