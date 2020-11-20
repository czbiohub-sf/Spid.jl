function filter_seqs_sites(seqs, filter_fun)
    filtered_sites = Vector{Vector{Char}}()
    observed_alleles = Set(Char[])
    for pos in 1:length(seqs[1])
        column = [seq[pos] for seq in seqs]
        if filter_fun(column)
            push!(filtered_sites, column)
        end
    end
    return filtered_sites
end

function filter_fasta_sites(in_fasta, filter_fun)
    names = String[]
    seqs = String[]
    read_fasta!(in_fasta, names, seqs);

    for seq in seqs
        @assert length(seq) == length(seqs[1])
    end

    filtered_sites = filter_seqs_sites(seqs, filter_fun)
    new_seqs = [Char[] for s in seqs]
    for site in filtered_sites
        for (i, allele) in enumerate(site)
            push!(new_seqs[i], allele)
        end
    end

    return [
        FASTA.Record(name, LongCharSeq(join(seq)))
        for (name, seq) in zip(names, new_seqs)
    ]
end

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

function is_core_site(column)
    for a in column
        if a=='N' || a=='-'
            return false
        end
    end
    return true
end

function get_core_only_seqs(in_fasta)
    return filter_fasta_sites(in_fasta, is_core_site)
end

export filter_fasta_sites
export variants_only_fasta
export get_variant_only_seqs
export get_core_only_seqs
