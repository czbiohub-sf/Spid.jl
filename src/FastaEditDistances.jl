module FastaEditDistances

using Distributed

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

"""
Returns a matrix whose upper triangle is the edit difference,
and lower triangle pairwise shared genome length

genome_arr sould be a 2d Character Array whose rows are samples,
columns are bases
"""
function compute_diffs_helper(genome_arr)
    n = size(genome_arr, 1)
    ret = zeros(Int, n, n)
    for l in 1:size(genome_arr, 2)
        skip_samples = [check_skip_allele(genome_arr[i, l]) for i in 1:n]
        for i in 1:(n-1)
            if !skip_samples[i]
                a_i = genome_arr[i, l]
                for j in (i+1):n
                    if !skip_samples[j]
                        a_j = genome_arr[j, l]
                        ret[i,j] += 1
                        if a_i != a_j
                            ret[j,i] += 1
                        end
                    end
                end
            end
        end
    end
    return ret
end

function parallel_compute_diffs(genome_arr, batch_size)
    L = size(genome_arr, 2)
    ret = @distributed (+) for i in range(1, stop=L, step=batch_size)
        compute_diffs_helper(genome_arr[:,i:min(L, i+batch_size-1)])
    end
    return ret
end

function seqs2arr(seqs)
    n = length(seqs)
    L = length(seqs[1])
    genome_arr = Array{Char}(undef, n, L)
    for i in 1:n
        for j in 1:L
            genome_arr[i,j] = seqs[i][j]
        end
    end
    return genome_arr
end

function fasta_edit_distances(fasta_name; batch_size=100)
    names = String[]
    seqs = String[]
    open(fasta_name) do fasta
        for record in FASTA.Reader(fasta)
            name = FASTA.identifier(record)
            seq = uppercase(String(sequence(CharSequence, record)))
            push!(names, name)
            push!(seqs, seq)
        end
    end

    res = parallel_compute_diffs(seqs2arr(seqs), batch_size)

    Contig1 = String[]
    Contig2 = String[]
    SharedGenomeLen = Int[]
    NumDiffs = Int[]

    for i = 1:(length(seqs)-1)
        for j = (i+1):length(seqs)
            push!(Contig1, names[i])
            push!(Contig2, names[j])
            push!(SharedGenomeLen, res[i,j])
            push!(NumDiffs, res[j,i])
        end
    end

    return DataFrame(
        Contig1 = Contig1,
        Contig2 = Contig2,
        SharedGenomeLen = SharedGenomeLen,
        NumDiffs = NumDiffs
    )
end

export fasta_edit_distances

end
