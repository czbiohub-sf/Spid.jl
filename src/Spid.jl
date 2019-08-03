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

include("FastaEditDistances.jl")
export fasta_edit_distances

include("VariantsOnlyFasta.jl")
export variants_only_fasta
export get_variant_only_seqs

include("Pileup2Consensus.jl")
export pileup_consensus_fasta
export bam2pileup2consensus

end # module
