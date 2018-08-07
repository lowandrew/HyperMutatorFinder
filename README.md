# HyperMutatorFinder

### Installation

Clone this repository, and install necessary python packages (you can `pip install -r requirements.txt`).

This is only intended for python 3.x, no python 2.x support.

Other depedencies (must be available on your $PATH):

- BBTools >= 38
- Samtools >= 1.6

### Usage

```
Attempts to find evolutionary rate of a strain, when provided with paired-end
FASTQ files and an assembly for that strain.

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        Full path to your FASTA-formatted assembly. Only works
                        with SPAdes assemblies right now. SKESA causes
                        breakpoints at most sites we're looking for. Other
                        assemblers have not been tested.
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        Full path to your forward (_R1) read file.
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        Full path to your reverse(_R2) read file.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Must not currently exist.
  -m MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        Minimum amount of coverage a site must have before it
                        will be considered for analysis. Defaults to 10.
  -w FILTER_DENSITY_WINDOW, --filter_density_window FILTER_DENSITY_WINDOW
                        Window size for identifying high-density SNV regions
                        that will be discarded from analysis. Defaults to 500,
                        meaning any region of 500 or fewer base pairs with
                        more than one SNV will not be considered.
  -c COVERAGE_FILTER, --coverage_filter COVERAGE_FILTER
                        Filter for coverage - sites that look heterogenous but
                        have higher than average coverage usually result from
                        two genes that are very close to identical that get
                        merged into one by SPAdes. This filter stops that.
                        Defaults to discarding if coverage is more than 1.5x
                        the average for the contig.
```

### Interpreting Results

WIP - so far, looks like anything that turns up with >30 multi-SNP positions is likely to be a HyperMutator, though
this may vary a lot by species.

