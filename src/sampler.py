import argparse
import random

def main():
    parser = argparse.ArgumentParser(
        description="Samples sequences from a FASTA file")
    parser.add_argument("-n",
                        type=int, 
                        help="the number of sequences to sample "
                        "(default: sample all)")
    parser.add_argument("-c",
                        type=int, 
                        help="the number of columns to sample "
                        "(default: sample all)")
    parser.add_argument("input_file",
                        type=argparse.FileType("r"),
                        help="the FASTA file to sample from")
    parser.add_argument("output_file",
                        type=argparse.FileType("w"),
                        help="the file to output the samples to")
    args = parser.parse_args()

    fasta_data = args.input_file.read().strip().split('\n')
    num_sequences = len(fasta_data) // 2
    seq_length = len(fasta_data[1])
    if args.n is None:
        args.n = num_sequences
    if args.c is None:
        args.c = seq_length
    if args.n <= 0:
        raise ValueError("Must sample at least 1 sequence")
    if args.n > num_sequences:
        raise ValueError(f"Cannot sample {args.n} sequences from a file with "
                         f"{num_sequences} sequences")
    if args.c <= 0:
        raise ValueError("Must sample at least 1 column")
    if args.c > seq_length:
        raise ValueError(f"Cannot sample {args.c} columns from sequences with "
                         f"length {seq_length}")
    seq_samples = random.sample(range(num_sequences), args.n)
    column_samples = sorted(random.sample(range(seq_length), args.c))
    for i in seq_samples:
        args.output_file.write(fasta_data[2*i] + '\n')
        seq = fasta_data[2*i+1]
        sampled_seq = ''.join(seq[j] for j in column_samples)
        args.output_file.write(sampled_seq + '\n')

if __name__ == "__main__":
    main()
