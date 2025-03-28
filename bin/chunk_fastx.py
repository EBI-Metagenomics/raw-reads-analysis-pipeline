import argparse
from collections import defaultdict
import gzip
import re
import math

parser = argparse.ArgumentParser(description='Chunk fasta and fastq files by base and read, handling single-end and paired-end reads. If both base and read targets are set then either of the criteria triggers the splitting of a chunk.')
parser.add_argument('-1', "--input_fp1", type=str,
                    required=True,
                    help="Input forward fasta/fastq filepath.")
parser.add_argument('-2', "--input_fp2", type=str,
                    required=False,
                    help="Input reverse fasta/fastq filepath.")
parser.add_argument('-o', "--output_pattern", type=str,
                    required=True,
                    help="Output pattern (will take directory, basename and extension and add chunk number and paired-end number).")
parser.add_argument('-b', "--base_count_target", type=str,
                    default='', required=False,
                    help="Target number of bases in each chunk (rounded down). K, M, G and T suffixes available (e.g. 1.5K)")
parser.add_argument('-r', "--read_count_target", type=str,
                    default='', required=False,
                    help="Target number of reads in each chunk (rounded down). K, M, G and T suffixes available (e.g. 1.5K).")
args = parser.parse_args()

def get_reads(infile, gz: bool):
    prev_line = None
    lines = []
    header = None
    base_count = None
    for line in infile:
        line_str = line.decode('utf8') if gz else line
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            if not header is None:
                yield (header, base_count, '\n'.join(lines))
            lines = []
            header, _ = re.findall(r"^(.*)([./](1|2))?$", line_str[1:].split()[0])[0][:2]
        if (not prev_line=='+') and (not line_str[0] in {'+','>','@'}):
            base_count = len(line_str.strip())
        lines.append(line_str.strip())
        prev_line = line_str.strip()
    else:
        yield (header, base_count, '\n'.join(lines))

def convert_size(size_str, size_suffix=None):
    if size_suffix is None:
        size_suffix = {v:i for i,v in enumerate(["", "K", "M", "G", "T"])}
    v, s = re.findall(r"^([0-9.,]+)([a-zA-Z]*)?$", size_str)[0]
    s = s[0] if s else ""
    return int(float(re.sub(r"[^0-9.]", '', v)) * math.pow(1000, size_suffix[s]))

if __name__ == '__main__':
    pe = bool(args.input_fp2)
    in_files = [args.input_fp1, args.input_fp2] if pe else [args.input_fp1]
    in_files = [(fp, fp[-3:]=='.gz') for fp in in_files]
    in_reads = [get_reads(gzip.open(fp, 'rb') if gz else open(fp, 'rt'), gz) for fp, gz in in_files]
    read_pairing_dict = defaultdict(set)
    read_stacks = [{} for _ in in_reads]
    base_counts = [0 for _ in in_reads]
    read_counts = [0 for _ in in_reads]
    stop_flag = [False for _ in in_reads]

    basename, extension = re.findall(r"^(.*)\.(f(ast)?[aq](\.gz)?)$", args.output_pattern)[0][:2]
    chunk_n = 1
    out_gz = extension[-3:]=='.gz'
    out_files = [None for _ in in_reads]

    size_suffix = {v:i for i,v in enumerate(["", "K", "M", "G", "T"])}
    read_count_target = convert_size(args.read_count_target, size_suffix) if args.read_count_target else False
    base_count_target = convert_size(args.base_count_target, size_suffix) if args.base_count_target else False

    while True:
        # read one read from each input file
        for i,g in enumerate(in_reads):
            if stop_flag[i]:
                continue
            try:
                k, base_count, v = next(g)
            except StopIteration:
                stop_flag[i] = True
                continue
            read_pairing_dict[k].add(i)
            read_stacks[i][k] = (base_count, v)

        # output reads that have been read in all input files
        add_k = [k for k,vs in read_pairing_dict.items() if len(vs)==len(in_reads)]
        for k in add_k:
            if any([v is None for v in out_files]):
                out_files = [f"{basename}{'.'+str(i+1) if pe else ''}.chunk-{chunk_n}.{extension}" for i,_ in enumerate(in_reads)]
                out_files = [gzip.open(fp, 'wb') if out_gz else open(fp, 'wt') for fp in out_files]

            for i, (s, o) in enumerate(zip(read_stacks, out_files)):
                base_count, v = s[k]
                out_str = v+'\n'
                o.write(out_str.encode() if out_gz else out_str)
                read_counts[i] += 1
                base_counts[i] += base_count
            # check chunking, reset counts and update out_files
            chunk = False
            if base_count_target:
                for c in base_counts:
                    if c>=base_count_target:
                        chunk = True
            if read_count_target:
                for c in read_counts:
                    if c>=read_count_target:
                        chunk = True
            if chunk:
                base_counts = [0 for _ in in_reads]
                read_counts = [0 for _ in in_reads]
                for f in out_files:
                    f.close()
                chunk_n += 1
                out_files = [None for _ in in_reads]

        # remove added reads from buffer
        for k in add_k:
            del read_pairing_dict[k]
            for i,d in enumerate(read_stacks):
                del read_stacks[i][k]

        # exit condition
        if all(stop_flag):
            break

    for f in out_files:
        if not f is None:
            f.close()

    print(f"Reads remaining: {[len(d) for d in read_stacks]}")

