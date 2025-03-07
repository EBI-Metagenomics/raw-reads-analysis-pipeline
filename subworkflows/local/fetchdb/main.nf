include { FETCHUNZIP } from "${projectDir}/modules/local/fetchunzip/main"

workflow FETCHDB {

    take:
        fetch_ch  // channel: [ val(meta), Map ]
        cache_path  // value channel
    main:
        // check if file is remote
        remote_ch = fetch_ch
            .filter{ meta,m -> m.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
            .map{ meta,m -> [meta,m,m.remote_path] }
        local_ch = fetch_ch
            .filter{ meta,m -> m.remote_path =~ /^[a-zA-Z]{2,}:\/\// }
            .filter{ meta,m -> m.local_path =~ /^[a-zA-Z]{2,}:\/\// }
            .map{ meta,m -> [meta,m.local_path] }
        // local_ch.view{ "local_ch - ${it}" }
        
        cache_path_ch = remote_ch.map{ 
            meta,m,fp -> 
            def fn = m.remote_path.tokenize('/').last()
            [meta, m, file("${cache_path}/${meta.id}/${fn}")] 
        }
        download_ch = cache_path_ch
            .filter{ meta,m,cache_fp -> !cache_fp.exists() }
            .map{ meta,m,cache_fp -> [meta,m] }
        cache_ch = cache_path_ch
            .filter{ meta,m,cache_fp -> cache_fp.exists() }
            .map{ meta,m,cache_fp -> [meta,cache_fp] } 
        // cache_ch.view{ "cache_ch - ${it}" }

        FETCHUNZIP(download_ch)
        downloaded_ch = FETCHUNZIP.out
        // downloaded_ch.view{ "downloaded_ch - ${it}" }
    emit:
        // join downloaded
        dbs = local_ch.mix(cache_ch).mix(downloaded_ch)
}
