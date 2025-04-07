import argparse
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

from alignment import Alignment

def main():
    parser = argparse.ArgumentParser(description="FastTree implemented in Python")
    parser.add_argument("input_file",
                        type=argparse.FileType('r'),
                        help="the aligned nucleotide sequences in fasta format")

    args = parser.parse_args()
    logger.info(f"Loading fasta file: {args.input_file.name}")
    fasta_data = args.input_file.read().strip().split('\n')

    alignment_dict = dict()
    for label_line, seq in zip(fasta_data[::2], fasta_data[1::2]):
        label = label_line[1:].split(' ', 1)[0]
        alignment_dict[label] = seq

    logger.info(f"Constructing profile matrices")
    alignment = Alignment(alignment_dict)
    logger.info(f"Profile matrices of {alignment.alignment_size} sequences "
        f"of length {alignment.alignment_length} successfully constructed")

if __name__ == "__main__":
    main()
