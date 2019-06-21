module FastaEditDistances

using Base.Threads
using BioSequences
using DataFrames

const valid_alleles = "ACGT";
const skip_chars = "N-";

function check_skip_allele(a)
    if occursin(a, skip_chars)
        return true
    elseif occursin(a, valid_alleles)
        return false
    else
        throw(DomainError(a, string("$a not a valid allele")))
    end
end

function count_seq_diffs(seq1, seq2)
    if length(seq1) != length(seq2)
        throw(ErrorException(string(
            "contig lengths differ: ",
            "$(length(seq1)) != $(length(seq2))")))
    end

    total_bases = 0
    n_diffs = 0
    n_skipped = 0

    seq1 = uppercase(seq1)
    seq2 = uppercase(seq2)
    for (a1, a2) in zip(seq1, seq2)
        total_bases += 1
        if check_skip_allele(a1) || check_skip_allele(a2)
            n_skipped += 1
        elseif a1 != a2
            n_diffs += 1
        end
    end
    return [total_bases-n_skipped, n_diffs]
end

function fasta_edit_distances(fasta_name)
    names = []
    seqs = []
    open(fasta_name) do fasta
        for record in FASTA.Reader(fasta)
            name = FASTA.identifier(record)
            seq = String(sequence(CharSequence, record))
            push!(names, name)
            push!(seqs, seq)
        end
    end

    jobs = []
    for i = 1:(length(seqs)-1)
        for j = (i+1):length(seqs)
            push!(jobs, [i j])
        end
    end

    results = Array{Integer}(undef, length(jobs), 2)
    Threads.@threads for idx = 1:length(jobs)
        i = jobs[idx][1]
        j = jobs[idx][2]
        results[idx, :] = count_seq_diffs(seqs[i], seqs[j])
    end

    return DataFrame(
        Contig1 = [names[j[1]] for j in jobs],
        Contig2 = [names[j[2]] for j in jobs],
        SharedGenomeLen = results[:, 1],
        NumDiffs = results[:, 2]
    )
end

export fasta_edit_distances

end
