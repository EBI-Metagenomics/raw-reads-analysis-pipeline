import argparse
import gzip
import typing

parser = argparse.ArgumentParser(description='Generates sample-sheet from ENA accession.')
parser.add_argument('-f', "--input_fp1", type=str,
                    required=True,
                    help="Input forward fasta/fastq filepath.")
parser.add_argument('-r', "--input_fp2", type=str,
                    required=True,
                    help="Input reverse fasta/fastq filepath.")
parser.add_argument('-o', "--output_fp1", type=str,
                    required=True,
                    help="Output forward fasta/fastq filepath.")
parser.add_argument('-l', "--output_fp2", type=str,
                    required=True,
                    help="Output reverse fasta/fastq filepath.")
parser.add_argument('-d', "--suffix_delimeter", type=str,
                    default='/', required=False,
                    help="Delimeter for suffix (default: /)")
parser.add_argument('-e', "--encoding", type=str,
                    default='utf8', required=False,
                    help="Encoding for gzipped files.")
parser.add_argument('-z', "--zip-output", dest='output_zip',
                    action='store_true', help="Output to gzipped file.")
parser.add_argument('-c', "--no-zip-output", dest='output_zip',
                    action='store_false', help="Don't output to gzipped file.")

args = parser.parse_args()

def rename_file_gzip(
    infile: gzip.GzipFile,
    outfile: gzip.GzipFile,
    suffix: str,
    encoding: str = 'utf8'
):
    prev_line = None
    for line in infile:
        line_str = line.decode(encoding)
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            newline = line_str.split()[0]
            newline += args.suffix_delimeter + suffix + '\n'
            outfile.write(newline.encode())
        else:
            outfile.write(line)
        prev_line = line_str.strip()

def rename_file_text(
    infile: typing.TextIO,
    outfile: typing.TextIO,
    suffix: str,
):
    prev_line = None
    for line in infile:
        if (not prev_line=='+') and (line[0] in {'>','@'}):
            newline = line.split()[0]
            newline += args.suffix_delimeter + suffix + '\n'
            outfile.write(newline)
        else:
            outfile.write(line)
        prev_line = line.strip()

if __name__ == '__main__':
    in_gzip_format = args.input_fp1[-3:] == '.gz'
    out_gzip_format = (args.output_fp1[-3:] == '.gz') if output_gzip is None else output_gzip

    # this should be replaced with better use of file objects
    if in_gzip_format:
        with gzip.open(args.input_fp1, 'rb') as in_f:
            if out_gzip_format:
                with gzip.open(args.output_fp1, 'wb') as out_f:
                    rename_file_gzip(in_f, out_f, '1', encoding=args.encoding)
            else:
                with open(args.output_fp1, 'wt') as out_f:
                    rename_file_text(in_f, out_f, '1')

        with gzip.open(args.input_fp2, 'rb') as in_f:
            if out_gzip_format:
                with gzip.open(args.output_fp2, 'wb') as out_f:
                    rename_file_gzip(in_f, out_f, '2', encoding=args.encoding)
            else:
                with open(args.output_fp2, 'wt') as out_f:
                    rename_file_text(in_f, out_f, '2')
    else:
        with open(args.input_fp1, 'rt') as in_f:
            if out_gzip_format:
                with gzip.open(args.output_fp1, 'wb') as out_f:
                    rename_file_gzip(in_f, out_f, '1', encoding=args.encoding)
            else:
                with open(args.output_fp1, 'wt') as out_f:
                    rename_file_text(in_f, out_f, '1')

        with open(args.input_fp2, 'rt') as in_f:
            if out_gzip_format:
                with gzip.open(args.output_fp2, 'wb') as out_f:
                    rename_file_gzip(in_f, out_f, '2', encoding=args.encoding)
            else:
                with open(args.output_fp2, 'wt') as out_f:
                    rename_file_text(in_f, out_f, '2')
