using Spid
using CSV
using Test
using BioSequences

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

rm("scratch", force=true, recursive=true)
mkdir("scratch")

merge_alignments(
    "scratch/test", ["Sample1.fa", "Sample2.fa", "Sample3.fa"],
    Regex("(?:.*/)?([^/]*)\\.fa")
)

for fname in ["test.fa", "test.core.fa",
              "test.variantsOnly.fa", "test.core.variantsOnly.fa",
              "test.fa.pairwise_diffs.csv", "test.core.fa.pairwise_diffs.csv"]
    open(fname) do f1
        open("scratch/$fname") do f2
            @test read(f1, String) == read(f2, String)
        end
    end
end

# test empty core genome doesn't crash
merge_alignments(
    "scratch/test", ["Sample1.fa", "Sample2.fa", "Sample3.fa", "SampleN.fa"],
    Regex("(?:.*/)?([^/]*)\\.fa")
)
