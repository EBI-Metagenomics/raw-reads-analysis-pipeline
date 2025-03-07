process ACCESSION2SAMPLESHEET {
    tag "${accession}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/python:3.11':
    //    'biocontainers/python:3.11' }"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    val(accession)

    output:
    path("${accession}_samplesheet.csv") 

    script:
    def out_fp = "${accession}_samplesheet.csv" 
    """
    #!/bin/bash
    python ${workflow.projectDir}/bin/generate_sample_sheet.py -a "${accession}" -r "read_run" -o "${out_fp}"
    """
}

// process JOINSAMPLESHEETS {
//     tag "${fps.join(', ')}"
//     label 'process_single'
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/python:3.11':
//         'biocontainers/python:3.11' }"
//     publishDir "${params.outdir}", mode: 'copy'
// 
//     input:
//     val(fps)
// 
//     output:
//     path(out_fp)
// 
//     script:
//     // concatenate files, skipping the header
//     def fps_str = fps.join(',')
//     def out_fp = 'samplesheet.csv'
//     """
//     #!/bin/python
//     with open("${out_fp}", 'wt') as out_f:
//         for fp in "${fps_str}".split(','):
//             with open(fp.strip(), 'rt') as in_f:
//                 out_f.writelines([l for l in in_f if len(l)>1][1:])
//     """
// }
