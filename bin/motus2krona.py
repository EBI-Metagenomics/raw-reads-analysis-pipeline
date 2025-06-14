import argparse
import fileinput
from collections import defaultdict

parser = argparse.ArgumentParser(description='Converts mOTUs output to krona tsv table. Needs -p and -q mOTUs flags to generate compatible output.')
parser.add_argument('-i', "--input_fp", type=str,
                    default='-', required=False,
                    help="Input mOTUs `.out` filepath. Use '-' for STDIN (default).")
parser.add_argument('-o', "--output_fp", type=str,
                    default='-', required=False,
                    help="Output TSV filepath. Use '-' for STDOUT (default).")
parser.add_argument('-m', "--motu", action='store_true',
                    help="Add mOTU as final taxonomic level (default: False)")
parser.add_argument('-t', "--tax_level", type=int,
                    default=None, required=False,
                    help="Taxonomic level (default: -1, i.e. all)")

args = parser.parse_args()

outfile = None if args.output_fp=='-' else open(args.output_fp, 'wt')

if __name__ == '__main__':
    parsed_data = defaultdict(float)
    for line in fileinput.input([] if args.input_fp=='-' else args.input_fp):
        if line[0] == '#':
            continue
        motu, tax, ncbi, count = line.strip().split('\t')
        count_v = float(count)
        if count_v==0.:
            continue
        tax_list = tax.split('|')[:args.tax_level]
        if args.motu:
            tax_list += [motu]
        parsed_data[tuple(tax_list)] += float(count_v)

    for tax_list, count in parsed_data.items():
        print('\t'.join((str(count),) + tax_list), file=outfile)



