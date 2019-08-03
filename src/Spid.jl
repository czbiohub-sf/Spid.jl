module Spid

using Base.Threads
using BioSequences
using CodecZlib
using DataFrames

const valid_alleles = "ACGT";
const skip_chars = "N-";

# Only used by script entry-points, but include them here to
# pre-compile
using ArgParse
using CSV

function read_fasta!(fasta, names, seqs)
    for record in FASTA.Reader(fasta)
        name = FASTA.identifier(record)
        seq = String(sequence(CharSequence, record))
        push!(names, name)
        push!(seqs, seq)
    end
end

function filter_fasta_sites(in_fasta, filter_fun)
    names = String[]
    seqs = String[]
    read_fasta!(in_fasta, names, seqs);

    for seq in seqs
        @assert length(seq) == length(seqs[1])
    end

    filtered_sites = filter_seqs_sites(seqs, filter_fun)
    new_seqs = [Char[] for s in seqs]
    for site in filtered_sites
        for (i, allele) in enumerate(site)
            push!(new_seqs[i], allele)
        end
    end

    return [
        FASTA.Record(name, CharSequence(join(seq)))
        for (name, seq) in zip(names, new_seqs)
    ]
end

function filter_seqs_sites(seqs, filter_fun)
    filtered_sites = Vector{Vector{Char}}()
    observed_alleles = Set(Char[])
    for pos in 1:length(seqs[1])
        column = [seq[pos] for seq in seqs]
        if filter_fun(column)
            push!(filtered_sites, column)
        end
    end
    return filtered_sites
end

export filter_fasta_sites

include("FastaEditDistances.jl")
export fasta_edit_distances

include("VariantsOnlyFasta.jl")
export variants_only_fasta
export get_variant_only_seqs

include("Pileup2Consensus.jl")
export pileup_consensus_fasta
export bam2pileup2consensus

end # module
