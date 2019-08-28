#!/usr/bin/env julia
using ArgParse
using Spid

s = ArgParseSettings()
@add_arg_table s begin
    "align_short_reads"
    help = "Align short reads with minimap2 & generate consensus fasta."
    action = :command

    "align_assembly"
    help = "Align assembly with minimap2 & generate re-aligned fasta."
    action = :command

    "merge_alignments"
    help = "Merge multiple alignments and compute summaries."
    action = :command
end

@add_arg_table s["align_short_reads"] begin
    "out_prefix"
    help = "Prefix for output files"
    arg_type = String
    required = true

    "ref"
    help = "Reference fasta. Must be bgzipped and indexed by samtools faidx."
    arg_type = String
    required = true

    "fastqs"
    help = "Fastq files to align to reference."
    nargs = '+'
    arg_type = String
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

    "--threads"
    help = "Threads. Note: currently only used in alignment stage, but not in the pileup stage."
    arg_type = Int
    default = 3
end

@add_arg_table s["align_assembly"] begin
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

@add_arg_table s["merge_alignments"] begin
    "out_prefix"
    help = "Prefix for output files"
    arg_type = String
    required = true

    "fastas"
    help = "Fasta files to merge."
    nargs = '+'
    arg_type = String
    required = true

    "--fasta2sample_regex"
    help = "Regex for fasta path. Sample name is the first capturing group."
    default = "(?:.*/)?([^/]*)\\.fa"
end


parsed_args = parse_args(s)
cmd = parsed_args["%COMMAND%"]

if cmd == "align_short_reads"
    parsed_args = parsed_args[cmd]

    align_short_reads(
        parsed_args["ref"], parsed_args["fastqs"],
        string(parsed_args["out_prefix"], ".bam"),
        string(parsed_args["out_prefix"], ".fa"),
        min_ac=parsed_args["min_ac"], min_af=parsed_args["min_af"],
        max_dp=parsed_args["max_dp"],
        threads=parsed_args["threads"]
    )
elseif cmd == "align_assembly"
    parsed_args = parsed_args[cmd]

    align_assembly(
        parsed_args["ref"], parsed_args["query"],
        string(parsed_args["out_prefix"], ".bam"),
        string(parsed_args["out_prefix"], ".fa"),
        preset=parsed_args["preset"], min_ac=parsed_args["min_ac"],
        min_af=parsed_args["min_af"], max_dp=parsed_args["max_dp"],
        threads=parsed_args["threads"]
    )
elseif cmd == "merge_alignments"
    parsed_args = parsed_args[cmd]

    merge_alignments(
        parsed_args["out_prefix"],
        parsed_args["fastas"],
        Regex(parsed_args["fasta2sample_regex"])
    )
else
    println("Unrecognized command $cmd")
end
