# SlowTree

A Python implementation of the FastTree algorithm for phylogenetic tree construction. See `writeup.pdf` for background and implementation details.

### Installation

Clone repository:

```shell
git clone https://github.com/dmiao623/SlowTree.git
```

Create a Python virtual environment (optional):

```shell
python -m venv venv/
source venv/bin/activate
```

Add required packages:

```
pip install numpy
pip install blosum
pip install newick
```

### Usage

First, obtain a multiple alignment of nucleotide or peptide sequences in FASTA file format. For example, [a multiple alignment of 16S ribosomal rRNA from the greengenes database](https://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta).

If desired, run the `sampler.py` script to obtain a smaller subset of this dataset for testing. For example,
```
python src/sampler.py -n 100 core_set_aligned.fasta sampled_core_set.fasta
```

Before running the algorithm, ensure that the correct constants file is imported in `constants.py`. `_constants.dna` should be imported when using a dataset of nucleotide sequences, and `_constants.peptide` should be imported when using a dataset of peptide sequences. 

To run the actual algorithm, type
```
python src/main.py --algo slowtree sampled_core_set.fasta tree.txt
```
This will run the SlowTree algorithm on the sequences in `sampled_core_set.fasta` and output the resulting tree in [Newick format](https://en.wikipedia.org/wiki/Newick_format) to `tree.txt`.

For more usage options, type
```
python src/main.py -h
```

### References

- Price M. N., Dehal P. S., & Arkin A. P. (2009).  
  **[FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix](https://academic.oup.com/mbe/article/26/7/1641/1128976)**. *Molecular Biology and Evolution, 26*(7), 1641–1650. https://doi.org/10.1093/molbev/msp077 

- Bruno W. J., Socci N. D., & Halpern A. L. (2000).  
  **[Weighted Neighbor Joining: A Likelihood‑Based Approach to Distance‑Based Phylogeny Reconstruction](https://academic.oup.com/mbe/article/17/1/189/975625)**. *Molecular Biology and Evolution, 17*(1), 189–197. https://doi.org/10.1093/oxfordjournals.molbev.a026231

