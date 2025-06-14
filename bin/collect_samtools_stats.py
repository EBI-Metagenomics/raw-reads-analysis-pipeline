from collections import defaultdict
import json
import argparse
parser = argparse.ArgumentParser(description='Generates sample-sheet from ENA accession.')
parser.add_argument('-i', "--input_files", nargs='+', default=[],
                    help="Input samtools stats files.")
parser.add_argument('-p', "--prefix", default='stats',
                    help="Prefix for output files.")

args = parser.parse_args()


if __name__ == '__main__':
    summary_results = defaultdict(lambda : defaultdict(dict))
    read_lengths = defaultdict(lambda : defaultdict(dict))
    mapping_qualities = defaultdict(lambda : defaultdict(dict))
    for fp in args.input_files:
        with open(fp, 'rt') as f:
            for l in f:
                if l[0] == '#':
                    continue
                l_split = l.strip().split('\t')
                if l_split[0] == 'SN':
                    _, k, v = l_split[:3]
                    summary_results[fp][k[:-1]] = float(v)
                if l_split[0] == 'RL':
                    _, k, v = l_split[:3]
                    read_lengths[fp][int(k)] = int(v)
                if l_split[0] == 'MAPQ':
                    _, k, v = l_split[:3]
                    mapping_qualities[fp][int(k)] = int(v)
    with open(f"{args.prefix}_summary_results.json", 'wt') as f:
        json.dump(summary_results, f)
    with open(f"{args.prefix}_read_lengths.json", 'wt') as f:
        json.dump(read_lengths, f)
    with open(f"{args.prefix}_mapping_qualities.json", 'wt') as f:
        json.dump(mapping_qualities, f)

