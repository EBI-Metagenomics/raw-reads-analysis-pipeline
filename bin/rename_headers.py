import argparse
import fileinput

parser = argparse.ArgumentParser(description='Generates sample-sheet from ENA accession.')
parser.add_argument('-i', "--input_fp", type=str,
                    default='-', required=False,
                    help="Input fasta/fastq filepath. Use '-' for STDIN (default).")
parser.add_argument('-o', "--output_fp", type=str,
                    default='-', required=False,
                    help="Output fasta/fastq filepath. Use '-' for STDOUT (default).")
parser.add_argument('-s', "--suffix", type=str,
                    default='', required=False,
                    help="Add suffixes (for paired-end reads) (default: False)")
parser.add_argument('-d', "--suffix_delimeter", type=str,
                    default='/', required=False,
                    help="Delimeter for suffix (default: /)")

args = parser.parse_args()

outfile = None if args.output_fp=='-' else open(args.output_fp, 'wt')

if __name__ == '__main__':
    prev_line = None
    for line in fileinput.input([] if args.input_fp=='-' else args.input_fp):
        if (not prev_line=='+') and (line[0] in {'>','@'}):
            newline = line.split()[0]
            if args.suffix:
                newline += args.suffix_delimeter + args.suffix
            print(newline.strip(), file=outfile)
        else:
            print(line.strip(), file=outfile)
        prev_line = line.strip()

