include { FETCHUNZIP } from '../../../modules/local/fetchunzip/main'

workflow FETCHDB {
    take:
    fetch_ch    // channel: val(meta)
    cache_path  // String

    main:
    remote_ch = fetch_ch
        .filter { meta -> meta.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
        .map { meta -> [meta, file(meta.remote_path)] }
    local_ch = fetch_ch
        .filter { meta -> meta.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
        .filter { meta -> meta.local_path =~ /^[a-zA-Z]{2,}:\/\// }
        .map { meta -> [meta, file(meta.local_path, checkIfExists: true)] }
    // local_ch.view { "local_ch - ${it}" }

    cache_path_ch = remote_ch.map { meta, _fp ->
        [meta, file("${projectDir}/${cache_path}/${meta.id}")]
    }
    // cache_path_ch.view { "cache_path_ch - ${it}" }
    download_ch = cache_path_ch
        .filter { _meta, cache_fp -> !cache_fp.exists() }
        .map { meta, _cache_fp -> [meta, meta.id, file(meta.remote_path, checkIfExists: true)] }
    cache_ch = cache_path_ch
        .filter { _meta, cache_fp -> cache_fp.exists() }
        .map { meta, cache_fp -> [meta, cache_fp] }
    // cache_ch.view { "cache_ch - ${it}" }
    download_ch.view { "download_ch - ${it}" }

    FETCHUNZIP(download_ch)
    downloaded_ch = FETCHUNZIP.out

    emit:
    dbs = local_ch.mix(cache_ch).mix(downloaded_ch)
}
