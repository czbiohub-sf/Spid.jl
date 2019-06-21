module Pileup2Consensus

#using ArgParse
using BioSequences
using CodecZlib

#function parse_commandline()
#    s = ArgParseSettings()
#
#    @add_arg_table s begin
#        "bam"
#            help = "Bam file. Must be sorted and indexed."
#            required = true
#        "ref"
#            help = "Reference fasta. Must be bgzipped and indexed by samtools faidx."
#            required = true
#        "--min_ac"
#            help = "Minimum read count of consensus allele."
#            arg_type = Int
#            default = 10
#        "--min_af"
#            help = "Minimum frequency of consensus allele."
#            arg_type = Float64
#            default = 0.9
#        "--max_dp"
#            help = "Maximum read depth; positions with higher depth are masked."
#            arg_type = Int
#            default = typemax(Int)
#    end
#
#    return parse_args(s)
#end
#
#function main()
#    parsed_args = parse_commandline()
#
#    ref_filename = parsed_args["ref"]
#    bam_filename = parsed_args["bam"]
#
#    check_samtools_version()
#    cmd = `samtools mpileup -f $ref_filename $bam_filename`
#    open(cmd) do pileup_stream
#        open(ref_filename) do ref_stream
#            ref_fasta = FASTA.Reader(GzipDecompressorStream(ref_stream))
#
#            consensus_fasta = pileup_consensus_fasta(
#                pileup_stream, ref_fasta,
#                parsed_args["min_ac"], parsed_args["min_af"], parsed_args["max_dp"]
#            )
#
#            w = FASTA.Writer(stdout)
#            for record in consensus_fasta
#                write(w, record)
#            end
#        end
#    end
#end

function check_samtools_version()
    open(`samtools --version`) do cmd
        version = match(r"samtools (\d+\.\d+)", readline(cmd)).captures[1]
        version = parse(Float64, version)
        min_version = 1.4
        recommended_version = 1.9
        if version < min_version
            throw(ErrorException(string(
                "samtools version $version too low.",
                " samtools>=$min_version required;",
                " samtools>=$recommended_version recommended.")))
        elseif version < recommended_version
            @warn "samtools=$version. samtools>=$recommended_version recommended."
        end
    end
end

function pileup_consensus_fasta(pileup_stream, ref_fasta,
                                ac_lower, af_lower, dp_upper)
    consensus_seq_list = initialize_consensus_seq(ref_fasta)
    consensus_seq_dict = Dict{String, Array{Char}}(
        name => seq for (name, seq) in consensus_seq_list)

    chrom, pos = "", 0

    allele_array = Char[]
    allele_counter = Dict{Char, Int}(a => 0 for a in "ACGTN-")
    for line in eachline(pileup_stream)
        chrom, pos, ref, dp, read_bases, base_quals = parse_pileup_line(
            line, chrom, pos)

        consensus_seq_dict[chrom][pos] = consensus_allele!(
            allele_counter, allele_array,
            read_bases, ref,
            ac_lower, af_lower, dp_upper
        )
    end

    return [
        FASTA.Record(name, CharSequence(join(seq)))
        for (name, seq) in consensus_seq_list
    ]
end

function parse_pileup_line(line, prev_chrom, prev_pos)
    fields = split(line)

    chrom = fields[1]
    pos = parse(Int, fields[2])
    @assert chrom != prev_chrom || pos > prev_pos

    ref_idx = 3
    @assert length(fields[ref_idx]) == 1
    ref = fields[ref_idx][1]

    dp = parse(Int, fields[4])
    read_bases = fields[5]
    base_quals = fields[6]

    return (chrom, pos, ref, dp, read_bases, base_quals)
end

function initialize_consensus_seq(ref_fasta)
    consensus_seq = []
    for record in ref_fasta
        name = FASTA.identifier(record)
        seq = sequence(CharSequence, record)
        push!(consensus_seq, (name, ['N' for i in 1:length(seq)]))
    end
    return consensus_seq
end

function consensus_allele!(allele_counter, allele_array,
                           read_bases, ref, ac_lower, af_lower, dp_upper)
    fill_alleles_array!(allele_array, read_bases, ref)
    dp = count_alleles!(allele_counter, allele_array)
    major_allele, max_ac = max_allele_and_count(allele_counter)
    max_af = max_ac / dp

    if max_ac >= ac_lower && max_af >= af_lower && dp <= dp_upper
        return major_allele
    else
        return 'N'
    end
end

function ref_and_bases(pileup_line)
    return ref, read_bases
end

function count_alleles!(allele_counter, allele_array)
    for k in keys(allele_counter)
        allele_counter[k] = 0
    end
    dp = 0
    for a in allele_array
        allele_counter[a] += 1
        dp += 1
    end
    return dp
end

function max_allele_and_count(allele_counter)
    major_allele = 'N'
    max_ac = 0
    for (a, ac) in allele_counter
        if ac > max_ac
            major_allele = a
            max_ac = ac
        end
    end
    return major_allele, max_ac
end

function fill_alleles_array!(allele_array, read_bases, ref)
    resize!(allele_array, 0)
    idx = 1
    while idx <= length(read_bases)
        idx = process_next_symbol!(allele_array, read_bases, idx, ref)
    end
end

function process_next_symbol!(allele_array, read_bases, idx, ref)
    next_char = read_bases[idx]
    if next_char == '+' || next_char == '-'
        return indel_bases_end(read_bases, idx) + 1
    elseif next_char == '^'
        return idx + 2
    elseif next_char == '$' || next_char == '<' || next_char == '>'
        return idx + 1
    else
        push!(allele_array, get_allele(next_char, ref))
        return idx + 1
    end
end

function indel_bases_end(read_bases, indel_start)
    digits_str = match(r"\d+", read_bases, indel_start).match
    digits_num = parse(Int, digits_str)
    return indel_start + length(digits_str) + digits_num
end

function get_allele(symbol, ref)
    if symbol == '.' || symbol == ','
        return ref
    elseif symbol == '*'
        return '-'
    else
        return uppercase(symbol)
    end
end

export pileup_consensus_fasta

#main()

end
