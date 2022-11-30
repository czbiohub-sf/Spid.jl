# SPID: SNP Pipeline for Infectious Disease

This package implements a simple alignment and SNP-calling pipeline
for detecting viral and bacterial infectious disease outbreak
clusters.

## Installation

The easiest way to use spid is through the docker container `czbiohub/spid`.
Currently, this is the best way to run the tool so please use this option for now.


## Usage

Spid can be run from julia by importing it with `using Spid`, or from
the command-line as `spid.jl [command]`.

There are 3 main commands:
- `align_short_reads`: Align short reads against a reference.
- `align_assembly`: Align already assembled fasta against a reference.
- `merge_alignments`: Combine samples generated from
  `align_short_reads` and `align_assembly` into a single dataset, and
  compute pairwise distances between samples.

To see the help for these functions, run
`spid.jl [command] --help` from the command line, or
`?{command}` within julia.

### Running through docker

First pull the image from dockerhub:
```{sh}
docker pull czbiohub/spid:<version>
```

To run the command-line utility, with the current directory mounted
at `/data` inside the docker image, do:
```{sh}
docker run -v $PWD:/data czbiohub/spid spid.jl [...]
```

To run julia through the docker container, with the current directory mounted at `/data`, do:
```{sh}
docker run -it -v $PWD:/data czbiohub/spid julia
```


### Manual installation

If you prefer not to use the docker container, you can manually
install as follows:

1. Install dependencies:
   - minimap2 >= 2.15
   - samtools >= 1.9
   - julia >= 1.0
2. Add the top-level folder as a local julia package
   - You can do this by typing `]` and then `add .` in julia
3. (Optional) Add the executable `bin/spid.jl` to your $PATH.

## Implementation Details

The pipeline aligns samples against a reference genome using minimap2,
then uses samtools to perform an mpileup, and then runs julia code to
call a consensus allele at each position according to user-specified
filters. It also runs julia code to compute the SNP distances between
every pair of samples.

The pipeline is meant to work well in an online setting where samples
are added incrementally, without affecting previous results. Thus,
alleles are called individually rather than jointly across samples.
It also reports a whole-genome SNP distance between every pair, which
does not change as new samples are added, in addition to the more
traditional core-genome SNP distance which does change as samples are
added (since the core genome size shrinks with sample size).

This pipeline is only meant for viral and bacterial haploid
sequences. It will not work for diploid or polyploid eukaryotic
sequences, though it can also work for haploid fungal sequences
(e.g., during infection when they are in the spore stage).
