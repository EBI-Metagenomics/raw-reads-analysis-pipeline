import re
import gzip
from collections import defaultdict
import argparse
from typing import Optional

parser = argparse.ArgumentParser(description='Standardise FASTA and FASTQ files by keeping only the first item in the header, adding a paired-end suffix to each header and splitting paired-end reads into two files.')
parser.add_argument('-1', "--input_fp1", type=str,
                    required=True,
                    help="Input forward fasta/fastq filepath.")
parser.add_argument('-2', "--input_fp2", type=str,
                    required=False,
                    help="Input reverse fasta/fastq filepath.")
parser.add_argument('-o', "--output_pattern", type=str,
                    required=True,
                    help="Output pattern (will take directory, basename and extension and add chunk number and paired-end number).")
parser.add_argument('-u', "--output_unpaired", action='store_true',
                    help="Output a fastx file containing all reads that weren't paired up.")
parser.add_argument('-z', "--output_gzip", action='store_true',
                    help="Gzip the output. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-c', "--no_output_gzip", action='store_false', dest='output_gzip',
                    help="Don't gzip the output. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-p', "--paired_end", action='store_true',
                    help="Flag specifying if input (and therefore outputs) are paired-end. If so then two files with be creates with suffixes '_1' and '_2'.")
parser.add_argument('-b', "--io_batch", type=int,
                    default=1000,
                    help="Number of FASTX records in a read/write batch.")
args = parser.parse_args()


def get_reads(
        infile,
        gz: bool,
        pe: Optional[int] = None
):
    prev_line = None
    lines = []
    header = None
    pe_n = None
    for line in infile:
        line_str = line.decode('utf8') if gz else line
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            if not header is None:
                yield (header, pe_n, '\n'.join(lines))
            lines = []
            header, pe_n = re.findall(r"^(.*?)([./][12])?$", line_str[1:].split()[0])[0][:2]
            if (args.paired_end) and (len(pe_n)>0):
                pe_n = int(pe_n[1:])-1
            else:
                if pe is not None:
                    pe_n = pe
                else:
                    pe_n = None
            if pe_n is None:
                lines.append(f"{line_str[0]}{header}")
            else:
                lines.append(f"{line_str[0]}{header}/{pe_n+1}")
        else:
            lines.append(line_str.strip())
        prev_line = line_str.strip()
    else:
        yield (header, pe_n, '\n'.join(lines))

if __name__ == '__main__':
    pe = bool(args.paired_end) or (args.input_fp2 is not None)
    in_files = [fp for fp in [args.input_fp1, args.input_fp2] if fp is not None]
    in_files = [(fp, fp[-3:]=='.gz') for fp in in_files]
    in_reads = [get_reads(gzip.open(fp, 'rb') if gz else open(fp, 'rt'), gz, i) for i, (fp, gz) in enumerate(in_files)]
    read_pairing_dict = defaultdict(set)
    read_stacks = [{} for _ in range(2 if pe else 1)]
    stop_flag = [False for _ in in_reads]

    basename, extension = re.findall(r"^(.*?)\.(f(ast)?[aq](\.gz)?)$", args.output_pattern)[0][:2]
    out_gz = (args.output_gzip) if (args.output_gzip is not None) else (extension[-3:]=='.gz')
    if pe:
        out_files = [f"{basename}.{str(i+1)}.{extension}" for i in range(2)]
    else:
        out_files = [f"{basename}.{extension}"]
    out_files = [gzip.open(fp, 'wb') if out_gz else open(fp, 'wt') for fp in out_files]

    while True:
        # read one read from each input file
        for i,g in enumerate(in_reads):
            if stop_flag[i]:
                continue
            for _ in range(args.io_batch):
                try:
                    k, pe_n, v = next(g)
                    read_pairing_dict[k].add(pe_n)
                    read_stacks[pe_n][k] = v
                except StopIteration:
                    stop_flag[i] = True
                    continue

        # output reads that have been read in all input files
        add_k = [k for k,vs in read_pairing_dict.items() if len(vs)==(2 if pe else 1)]
        for i, (s, o) in enumerate(zip(read_stacks, out_files)):
            for k in add_k:
                v = s[k]
                out_str = v+'\n'
                o.write(out_str.encode() if out_gz else out_str)

        # remove added reads from buffer
        for k in add_k:
            del read_pairing_dict[k]
            for i,_ in enumerate(read_stacks):
                del read_stacks[i][k]

        # exit condition
        if all(stop_flag):
            break

    for f in out_files:
        if not f is None:
            f.close()

    if args.output_unpaired:
        for i, d in enumerate(read_stacks):
            out_file = f"{basename}.unpaired{'' if i is None else '_'+str(i)}.{extension}"
            out_file = gzip.open(out_file, 'wb') if out_gz else open(out_file, 'wt')
            for k, v in d.items():
                out_str = v+'\n'
                out_file.write(out_str.dencode() if out_gz else out_str)
