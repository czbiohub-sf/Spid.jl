using Spid
using CSV
using Test
using BioSequences

@test fasta_edit_distances("test.fa", batch_size=3) == CSV.read("test.dist.csv")

open("test.fa") do test_fa
    open("test.variants_only.fa") do test_variants_fa
        @test get_variant_only_seqs(test_fa) == [
            rec for rec in FASTA.Reader(test_variants_fa)]
    end
end

open("test.pileup") do pileup_f
    open("test.pileup.ref.fa") do ref_f
        open("test.pileup.out.fa") do out_f
            ref_fasta = FASTA.Reader(ref_f)
            @test pileup_consensus_fasta(
                pileup_f, ref_fasta, 10, .9, typemax(Int)) == [
                    rec for rec in FASTA.Reader(out_f)]
        end
    end
end
