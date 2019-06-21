using Spid
using CSV
using Test

@test fasta_edit_distances("test.fa") == CSV.read("test.dist.csv")
