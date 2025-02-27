import requests
import argparse
import pandas as pd
import os
import io

parser = argparse.ArgumentParser(description='Generates sample-sheet from ENA accession.')
parser.add_argument('-a', "--accession", type=str, help="ENA accession code")
parser.add_argument('-r', "--result", type=str, default='read_run', help="Type of result (read_run|sample|experiment)")
parser.add_argument('-o', "--out_fp", type=str, help="Output filepath for sample-sheet")
args = parser.parse_args()

base_url = 'https://www.ebi.ac.uk/ena/portal/api/filereport'# 'PRJEB51728'

def parse_reads_files(s:str):
    fps = [f'ftp://{v}' for v in s.split(';')]
    if len(fps)==1:
        return [fps[0],None,None]
    if len(fps)==2:
        return list(sorted(fps)) + [None]
    if len(fps)==3:
        return list(sorted(fps, key=lambda x:(-len(x),x)))
    if len(fps)>3:
        raise Exception(f"Number of files >3 ({len(fps)})\n{'\n'.join(fps)}")

def parse_md5s(s:str):
    md5s = list(s.split(';'))
    return [md5s[i] if i<len(md5s) else None for i in range(3)]

if __name__ == '__main__':
    resp = requests.get(base_url, params={'accession': args.accession, 'result': args.result, 'fields': ','.join(['run_accession','fastq_ftp','fastq_md5','fastq_bytes'])})
    if not resp.status_code==200:
        raise Exception(f'Repsonse code {resp.code}')
    r = pd.read_csv(io.StringIO(resp.text), sep='\t')

    r[['fastq_1', 'fastq_2', 'fastq_barcode']] = [parse_reads_files(v) for v in r.fastq_ftp]
    r[['fastq_1_md5', 'fastq_2_md5', 'fastq_barcode_md5']] = [parse_md5s(v) for v in r.fastq_md5]
    r = r.rename(columns={'run_accession': 'sample'})

    os.makedirs(os.path.dirname(os.path.abspath(args.out_fp)), exist_ok=True)
    r[['sample', 'fastq_1', 'fastq_2', 'fastq_barcode', 'fastq_1_md5', 'fastq_2_md5', 'fastq_barcode_md5']].to_csv(args.out_fp, index=False)

