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

include("FastaEditDistances.jl")
include("FilterFastaSites.jl")
include("Pileup2Consensus.jl")

end # module
