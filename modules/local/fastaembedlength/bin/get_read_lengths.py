import argparse
import gzip

parser = argparse.ArgumentParser(description='Add the length of the sequence to the header of the sequence in a fasta or fastq file.')
parser.add_argument('-i', "--input_fp", type=str,
                    required=True,
                    help="Input fasta/fastq filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    default='',
                    help="Filepath for output tsv containing read lengths.")
args = parser.parse_args()

def get_reads(infile, gz: bool):
    prev_line = None
    lines = []
    header = None
    counting = False
    count = 0
    for line in infile:
        line_str = line.decode('utf8') if gz else line
        line_str_strip = line_str.strip()
        if line_str[0]=='+':
            counting = False
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            if not header is None:
                l = sum([len(v) for v in lines[1:]])
                yield (header, count, '\n'.join(lines[1:]))
            lines = []
            header = line_str.split()[0]
            counting = True
            count = 0
        else:
            if counting:
                count += len(line_str_strip)
        lines.append(line_str_strip)
        prev_line = line_str_strip
    else:
        yield (header, count, '\n'.join(lines[1:]))

if __name__ == '__main__':
    gz = args.input_fp[-3:]=='.gz'
    in_file = gzip.open(args.input_fp, 'rb') if gz else open(args.input_fp, 'rt')
    out_file = open(args.output_fp, 'wt')

    for k, c, v in get_reads(in_file, gz):
        out_str = f"{k}\t{c}\n"
        out_file.write(out_str)

    in_file.close()
    out_file.close()

