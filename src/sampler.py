import argparse
import random

def main():
    parser = argparse.ArgumentParser(
        description="Samples sequences from a FASTA file")
    parser.add_argument("-n",
                        type=int, 
                        help="the number of sequences to sample",
                        required=True)
    parser.add_argument("input_file",
                        type=argparse.FileType("r"),
                        help="the FASTA file to sample from")
    parser.add_argument("output_file",
                        type=argparse.FileType("w"),
                        help="the file to output the samples to")
    args = parser.parse_args()

    if args.n <= 0:
        raise ValueError("Must sample at least 1 sequence")
    fasta_data = args.input_file.read().strip().split('\n')
    num_sequences = len(fasta_data) // 2
    if args.n > num_sequences:
        raise ValueError(f"Cannot sample {args.n} sequences from a file with "
                         f"{num_sequences} sequences")
    samples = random.sample(range(num_sequences), args.n)
    for i in samples:
        args.output_file.write(fasta_data[2*i] + '\n')
        args.output_file.write(fasta_data[2*i+1] + '\n')

if __name__ == "__main__":
    main()
