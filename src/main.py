import argparse
import logging
import os
import resource
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

from alignment import Alignment
from tree_builder import TreeBuilder
from benchmarks.neighbor_joining import neighbor_joining
from benchmarks.random_joining import random_joining
from math import isqrt
import newick
import time

def get_peak_mem_mb():
    # ru_maxrss is in kilobytes on Linux, bytes on macOS
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if os.uname().sysname == 'Darwin':
        return peak / (1024**2)
    else:
        return peak / 1024

def main():
    sys.setrecursionlimit(10_000)
    print("Recursion limit:", sys.getrecursionlimit())

    parser = argparse.ArgumentParser(description="FastTree implemented in Python")
    parser.add_argument("--algo",
                        type=str,
                        help="the algorithm used to construct the tree",
                        required=True,
                        choices=["nj", "random", "slowtree"])
    parser.add_argument("input_file",
                        type=argparse.FileType("r"),
                        help="the aligned nucleotide sequences in fasta format")
    parser.add_argument("output_file",
                        type=argparse.FileType("w"),
                        help="the file to output the tree to")

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

    time_elapsed = time.perf_counter()
    if args.algo == "nj":
        newick.dump(neighbor_joining(alignment), args.output_file)
    elif args.algo == "random":
        newick.dump(random_joining(alignment), args.output_file)
    else:
        tree_builder = TreeBuilder(alignment,
                                   refresh_interval=isqrt(alignment.alignment_size))
        newick.dump(tree_builder.build(), args.output_file)

    time_elapsed = time.perf_counter() - time_elapsed
    logger.info(f"Elapsed time: {time_elapsed:.3f} s")
    logger.info(f"Peak memory usage: {get_peak_mem_mb():.2f} MiB")

if __name__ == "__main__":
    main()
