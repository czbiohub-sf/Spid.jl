using Printf

println("Time to load Spid")
@time using Spid

n = 20
L = 1000000
small_arr = Base.rand("ACGT", (2, 2))
genome_arr = Base.rand("ACGT", (n, L));

println("Time to compile parallel_compute_diffs")
@time Spid.FastaEditDistances.parallel_compute_diffs(small_arr, 1)

@printf("Time to run parallel_compute_diffs with n=%d, L=%d\n", n, L)
@time Spid.FastaEditDistances.parallel_compute_diffs(genome_arr, 100)
