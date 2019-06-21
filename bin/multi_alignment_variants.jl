#!/usr/bin/env julia
using BioSequences
using Spid

w = FASTA.Writer(stdout)
for record in get_variant_only_seqs(stdin)
    write(w, record)
end
