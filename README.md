# Create mRNA FASTA from genome assembly and BED files
This command line interface package can be used to create exon or intron FASTA files for designing \
fluorescent in situ hybridization (FISH) probes or for other genomic purposes.

## Dependencies 
- python >= 3.7.0
- pandas >= 1.1.3
- click >= 7.1.2
- biopython >= 1.79

## Installation
Navaigate to fasta_file_creator directory, then run the following:

```
pip install .
```
## Function arguments
```
create_pw2_fasta --help

Options:
  -f, --fasta TEXT      FASTA File of organism cDNA
  -b, --bed TEXT        BED file containing genomic coordinates
  -gf, --gff TEXT       gff file of organism  [required]
  -bio, --biotype TEXT  what type of molecule is this? (ex.introns or exons)
                        [required]

  -s, --sequences TEXT  Probe sequence csv file
  -g, --genenames TEXT  Gene names csv file for probe sequences
  -fi, --filename TEXT  output file name  [required]
  --help                Show this message and exit.
```
