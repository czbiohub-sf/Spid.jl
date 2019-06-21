#!/usr/bin/env julia
using ArgParse
using BioSequences
using Spid

s = ArgParseSettings()

@add_arg_table s begin
    "bam"
        help = "Bam file. Must be sorted and indexed."
        required = true
    "ref"
        help = "Reference fasta. Must be bgzipped and indexed by samtools faidx."
        required = true
    "--min_ac"
        help = "Minimum read count of consensus allele."
        arg_type = Int
        default = 10
    "--min_af"
        help = "Minimum frequency of consensus allele."
        arg_type = Float64
        default = 0.9
    "--max_dp"
        help = "Maximum read depth; positions with higher depth are masked."
        arg_type = Int
        default = typemax(Int)
end

parsed_args = parse_args(s)

bam2pileup2consensus(
    parsed_args["bam"], stdout,
    parsed_args["ref"],
    parsed_args["min_ac"], parsed_args["min_af"], parsed_args["max_dp"]
)
