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

function merge_sample_fastas(samples, in_fastas, out_fasta)
    # TODO check that fastas have same contig names and lengths
    open(out_fasta, "w") do f
        w = FASTA.Writer(f)
        for (sample_name, sample_fasta) in zip(samples, in_fastas)
            write(w, FASTA.Record(sample_name,
                                  concat_fasta_contigs(sample_fasta)))
        end
    end
end

function merge_and_summarize_sample_fastas(
    out_prefix, fastas, samples)

    samples_dict = Dict{String, String}()
    if samples == nothing
        samples = fastas
    elseif typeof(samples) == Regex
        samples = [match(samples, fa)[1] for fa in fastas]
    end

    merged_fa = "$out_prefix.fa"
    merge_sample_fastas(samples, fastas, merged_fa)
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

export merge_and_summarize_sample_fastas

include("FastaEditDistances.jl")
include("FilterFastaSites.jl")
include("Pileup2Consensus.jl")
include("Align.jl")

end # module
