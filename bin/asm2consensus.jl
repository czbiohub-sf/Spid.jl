#!/usr/bin/env julia
using ArgParse
using Spid

s = ArgParseSettings()

@add_arg_table s begin
    "out_prefix"
        help = "Prefix for output files"
        arg_type = String
        required = true
    "ref"
    help = "Reference fasta. Must be bgzipped and indexed by samtools faidx."
    arg_type = String
    required = true

    "query"
    help = "Query assembly."
    arg_type = String
    required = true

    "--preset"
    help = "Preset for minimap2"
    arg_type = String
    default = "asm5"

    "--min_ac"
    help = "Minimum read count of consensus allele."
    arg_type = Int
    default = 1

    "--min_af"
    help = "Minimum frequency of consensus allele."
    arg_type = Float64
    default = 1.0

    "--max_dp"
    help = "Maximum read depth; positions with higher depth are masked."
    arg_type = Int
    default = 1

    "--threads"
    help = "Threads. Note: currently only used in alignment stage, but not in the pileup stage."
    arg_type = Int
    default = 3
end

parsed_args = parse_args(s)

asm2consensus(
    parsed_args["ref"], parsed_args["query"],
    string(parsed_args["out_prefix"], ".bam"),
    string(parsed_args["out_prefix"], ".fa"),
    preset=parsed_args["preset"], min_ac=parsed_args["min_ac"],
    min_af=parsed_args["min_af"], max_dp=parsed_args["max_dp"],
    threads=parsed_args["threads"]
)
