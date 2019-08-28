module Spid

using Base.Threads
using BioSequences
using CodecZlib
using DataFrames
using CSV

const valid_alleles = "ACGT";
const skip_chars = "N-";

# Only used by script entry-points, but include them here to
# pre-compile
using ArgParse

function read_fasta!(fasta, names, seqs)
    for record in FASTA.Reader(fasta)
        name = FASTA.identifier(record)
        seq = String(sequence(CharSequence, record))
        push!(names, name)
        push!(seqs, seq)
    end
end

function concat_fasta_contigs(fasta_path)
    open(fasta_path) do f
        char_seq = Char[]
        stream = f
        if endswith(fasta_path, ".gz")
            stream = GzipDecompressorStream(f)
        end
        for record in FASTA.Reader(stream)
            append!(char_seq, sequence(CharSequence, record))
        end
        return CharSequence(join(char_seq))
    end
end

function combine_sample_fastas(samples, in_fastas, out_fasta)
    # TODO check that fastas have same contig names and lengths
    open(out_fasta, "w") do f
        w = FASTA.Writer(f)
        for (sample_name, sample_fasta) in zip(samples, in_fastas)
            write(w, FASTA.Record(sample_name,
                                  concat_fasta_contigs(sample_fasta)))
        end
    end
end

"""
    merge_alignments(out_prefix, fastas, samples)

Merge several samples into a combined FASTA and compute distance matrix.

Outputs the following FASTA files:
- `{out_prefix}.fa`: all sites
- `{out_prefix}.core.fa`: core sites only (no missing alleles)
- `{out_prefix}.variantsOnly.fa`: variant sites only
- `{out_prefix}.core.variantsOnly.fa`: core variants only

Also outputs pairwise distances in the following CSV files:
- `{out_prefix}.fa.pairwise_diffs.csv`: all sites
- `{out_prefix}.core.fa.pairwise_diffs.csv`: core sites only

`fastas` should be a list of FASTA filenames aligned against the same
reference, e.g. via [`align_short_reads`](@ref) or
[`align_assembly`](@ref). They should have the same number of contigs,
with the same lengths, in the same order. The contig names are
ignored. The contigs for each sample are concatenated into a
single sequence in the output FASTA.

`samples` should be either a list of the sample names, or a `Regex`
for extracting a sample name from the FASTA filename; in particular,
the first capturing group of the regex is taken as the sample name.
"""
function merge_alignments(
    out_prefix, fastas, samples)

    samples_dict = Dict{String, String}()
    if samples == nothing
        samples = fastas
    elseif typeof(samples) == Regex
        samples = [match(samples, fa)[1] for fa in fastas]
    end

    merged_fa = "$out_prefix.fa"
    combine_sample_fastas(samples, fastas, merged_fa)
    CSV.write("$out_prefix.fa.pairwise_diffs.csv",
              fasta_edit_distances(merged_fa))

    variants_fa = "$out_prefix.variantsOnly.fa"
    open(merged_fa) do in_f
        open(variants_fa, "w") do out_f
            w = FASTA.Writer(out_f)
            for record in get_variant_only_seqs(in_f)
                write(w, record)
            end
        end
    end

    core_fa = "$out_prefix.core.fa"
    open(merged_fa) do in_f
        open(core_fa, "w") do out_f
            w = FASTA.Writer(out_f)
            for record in get_core_only_seqs(in_f)
                write(w, record)
            end
        end
    end
    CSV.write("$out_prefix.core.fa.pairwise_diffs.csv",
              fasta_edit_distances(core_fa))

    core_variants_fa = "$out_prefix.core.variantsOnly.fa"
    open(core_fa) do in_f
        open(core_variants_fa, "w") do out_f
            w = FASTA.Writer(out_f)
            for record in get_variant_only_seqs(in_f)
                write(w, record)
            end
        end
    end
end

export merge_alignments

include("FastaEditDistances.jl")
include("FilterFastaSites.jl")
include("Pileup2Consensus.jl")
include("Align.jl")

end # module
