function get_variant_only_seqs(in_fasta)
    return filter_fasta_sites(in_fasta, is_variant_site)
end

function is_variant_site(column)
    observed_alleles = Set(Char[])
    fill_observed_alleles!(observed_alleles, column);
    return (length(observed_alleles) >= 2)
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
