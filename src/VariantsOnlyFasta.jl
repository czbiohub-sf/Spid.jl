module VariantsOnlyFasta

using BioSequences

const valid_alleles = "ACGT";
const skip_chars = "N-";

function get_variant_only_seqs(in_fasta)
    names = String[]
    seqs = String[]
    read_fasta!(in_fasta, names, seqs);

    for seq in seqs
        @assert length(seq) == length(seqs[1])
    end

    variants = get_variants(seqs)
    new_seqs = [Char[] for s in seqs]
    for variant in variants
        for (i, allele) in enumerate(variant)
            push!(new_seqs[i], allele)
        end
    end

    return [
        FASTA.Record(name, CharSequence(join(seq)))
        for (name, seq) in zip(names, new_seqs)
    ]
end

function get_variants(seqs)
    variants = Vector{Vector{Char}}()
    observed_alleles = Set(Char[])
    for pos in 1:length(seqs[1])
        column = [seq[pos] for seq in seqs]
        empty!(observed_alleles)
        fill_observed_alleles!(observed_alleles, column);
        if length(observed_alleles) >= 2
            push!(variants, column)
        end
    end
    return variants
end

function read_fasta!(fasta, names, seqs)
    for record in FASTA.Reader(fasta)
        name = FASTA.identifier(record)
        seq = String(sequence(CharSequence, record))
        push!(names, name)
        push!(seqs, seq)
    end
end

function fill_observed_alleles!(observed_alleles, column)
    for a in column
        if occursin(a, valid_alleles)
            push!(observed_alleles, a);
        elseif !occursin(a, skip_chars)
            throw(DomainError(a, string("$a not a valid allele")));
        end
    end
end

function variants_only_fasta(in_fasta, out_fasta)
    variant_only_seqs = get_variant_only_seqs(in_fasta);
    w = FASTA.Writer(out_fasta)
    for record in variant_only_seqs
        write(w, record)
    end
end

export get_variant_only_seqs
export variants_only_fasta

end
