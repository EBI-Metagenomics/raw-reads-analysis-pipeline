import requests
import argparse
import numpy as np
import pandas as pd
import os
import io
import re

parser = argparse.ArgumentParser(description='Generates sample-sheet from ENA accession.')
parser.add_argument('-a', "--accession", type=str, help="ENA accession code")
parser.add_argument('-r', "--result", type=str, default='read_run', help="Type of result (read_run|sample|experiment)")
parser.add_argument('-o', "--out_fp", type=str, help="Output filepath for sample-sheet")
parser.add_argument("--s3", action='store_true', help="Convert paths to internal EMBL-EBI s3 paths")
parser.add_argument("--s3_bucket", type=str, default='era-public', help="Data bucket for s3 file protocol")

args = parser.parse_args()

base_url = 'https://www.ebi.ac.uk/ena/portal/api/filereport'
# fields = ['study_accession', 'run_accession', 'fastq_ftp', 'library_layout', 'library_strategy', 'platform']
fields = ['study_accession', 'run_accession', 'fastq_ftp', 'library_layout', 'library_strategy', 'instrument_platform']

def parse_reads_files(s : str):
    fps = [f'ftp://{v}' for v in s.split(';')]
    if len(fps)==1:
        return [fps[0],None,None]
    if len(fps)==2:
        return list(sorted(fps)) + [None]
    if len(fps)==3:
        return list(sorted(fps, key=lambda x:(-len(x),x)))
    if len(fps)>3:
        raise Exception(f"Number of files >3 ({len(fps)}); {', '.join(fps)}")

def parse_md5s(s : str):
    md5s = list(s.split(';'))
    return [md5s[i] if i<len(md5s) else None for i in range(3)]

def fetch_accession(accession : str):
    resp = requests.get(base_url, params={'accession': accession, 'result': args.result, 'fields': ','.join(fields)})
    if not resp.status_code==200:
        raise Exception(f'Response code {resp.status_code}')
    r = pd.read_csv(io.StringIO(resp.text), sep='\t')
    r[['library_layout', 'library_strategy', 'instrument_platform']] = np.array([r[k].map(lambda x: x.lower() if isinstance(x,str) else x) for k in ['library_layout', 'library_strategy', 'instrument_platform']]).T
    r[['fastq_1', 'fastq_2', 'fastq_barcode']] = [parse_reads_files(v) for v in r.fastq_ftp]
    # r[['fastq_1_md5', 'fastq_2_md5', 'fastq_barcode_md5']] = [parse_md5s(v) for v in r.fastq_md5]
    r = r.rename(columns={'run_accession': 'reads_accession', 'instrument_platform': 'platform'})
    r['fetch_accession'] = accession

    return r

def convert_fp(fp : str, protocol : str, bucket : str):
    prefix_pattern = r'ftp://ftp\.sra\.ebi\.ac\.uk/vol\d+/'
    fp_ = re.split(prefix_pattern, fp)[1]
    return f"{protocol}://{bucket}/{fp_}"

if __name__ == '__main__':
    accessions = [v.strip() for v in re.split('[,; ]', args.accession)]
    dfs = [fetch_accession(a) for a in accessions if a]
    r = pd.concat(dfs, axis=0)

    if args.s3:
        r[['fastq_1', 'fastq_2']] = np.array([list(r[k].map(lambda x: convert_fp(x, protocol='s3', bucket=args.s3_bucket)).values) for k in ['fastq_1', 'fastq_2']]).T

    os.makedirs(os.path.dirname(os.path.abspath(args.out_fp)), exist_ok=True)
    # r[['study_accession', 'run_accession', 'fastq_1', 'fastq_2', 'library_layout', 'library_strategy', 'platform']].to_csv(args.out_fp, index=False)
    r[['study_accession', 'reads_accession', 'fastq_1', 'fastq_2', 'library_layout', 'library_strategy', 'platform']].to_csv(args.out_fp, index=False)
