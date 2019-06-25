module Spid

# Only used by script entry-points, but include them here to
# pre-compile
using ArgParse
using CSV

include("FastaEditDistances.jl")
using .FastaEditDistances
export fasta_edit_distances

include("VariantsOnlyFasta.jl")
using .VariantsOnlyFasta
export variants_only_fasta
export get_variant_only_seqs

include("Pileup2Consensus.jl")
using .Pileup2Consensus
export pileup_consensus_fasta
export bam2pileup2consensus

end # module
