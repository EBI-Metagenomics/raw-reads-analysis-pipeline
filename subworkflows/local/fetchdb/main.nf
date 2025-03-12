include { FETCHUNZIP } from '../../../modules/local/fetchunzip/main'

workflow FETCHDB {
    take:
    fetch_ch // channel: [ val(meta), Map ]
    cache_path // value channel

    main:
    // check if file is remote
    remote_ch = fetch_ch
        .filter { _meta, db_map -> db_map.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
        .map { meta, db_map -> [meta, db_map, file(db_map.remote_path)] }
    local_ch = fetch_ch
        .filter { _meta, db_map -> db_map.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
        .filter { _meta, db_map -> db_map.local_path =~ /^[a-zA-Z]{2,}:\/\// }
        .map { meta, db_map -> [meta, file(db_map.local_path)] }
    // local_ch.view{ "local_ch - ${it}" }

    cache_path_ch = remote_ch.map { meta, db_map, _fp ->
        [meta, db_map, file("${cache_path}/${meta.id}")]
    }
    download_ch = cache_path_ch
        .filter { _meta, _db_map, cache_fp -> !cache_fp.exists() }
        .map { meta, db_map, _cache_fp -> [meta, db_map, file(db_map.remote_path)] }
    cache_ch = cache_path_ch
        .filter { _meta, _db_map, cache_fp -> cache_fp.exists() }
        .map { meta, _db_map, cache_fp -> [meta, cache_fp] }
    // cache_ch.view{ "cache_ch - ${it}" }

    FETCHUNZIP(download_ch)
    downloaded_ch = FETCHUNZIP.out

    emit:
    dbs = local_ch.mix(cache_ch).mix(downloaded_ch)
}
