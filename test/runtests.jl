using Spid
using CSV
using Test
using BioSequences

@test fasta_edit_distances("test.fa") == CSV.read("test.dist.csv")

open("test.fa") do test_fa
    open("test.variants_only.fa") do test_variants_fa
        @test get_variant_only_seqs(test_fa) == [
            rec for rec in FASTA.Reader(test_variants_fa)]
    end
end
