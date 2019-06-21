#!/usr/bin/env julia
using CSV
using Spid

in_fasta, out_csv = ARGS
CSV.write(out_csv, fasta_edit_distances(in_fasta))
