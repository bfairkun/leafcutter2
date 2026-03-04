import sys
import gzip
import argparse


def star2junc(input_file, output):
    """Convert STAR SJ.out.tab file to BED6 junction format.

    Parameters
    ----------
    input_file : str
        Path to STAR SJ.out.tab file.
    output : str
        Path to output gzipped BED file.
    """
    with open(input_file, 'r') as fh:
        with gzip.open(output, 'wb') as out_fh:
            for line in fh:
                row = line.rstrip().split('\t')
                chrom = row[0]
                start = str(int(row[1])-1)  # Coordinates are 1-based
                end = row[2]
                counts = row[6]
                if row[3] == '1':
                    strand = '+'
                else:
                    strand = '-'
                out_line = ('\t'.join([chrom, start, end, '.', counts, strand]) + '\n').encode()
                out_fh.write(out_line)


def main_cli():
    parser = argparse.ArgumentParser(
        description="Convert STAR SJ.out.tab file to BED6 junction format."
    )
    parser.add_argument("input_file", help="STAR SJ.out.tab file")
    parser.add_argument("output", help="Output gzipped BED file")
    args = parser.parse_args()
    star2junc(args.input_file, args.output)


if __name__ == "__main__":
    main_cli()
