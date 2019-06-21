module Spid

include("FastaEditDistances.jl")
using .FastaEditDistances
export fasta_edit_distances

include("VariantsOnlyFasta.jl")
using .VariantsOnlyFasta
export variants_only_fasta
export get_variant_only_seqs

end # module
