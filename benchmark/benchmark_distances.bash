#/bin/bash
set -euxo pipefail
{
    export JULIA_PROJECT=$(realpath ..)
    julia benchmark_distances.jl
    julia -p 1 benchmark_distances.jl
    julia -p 2 benchmark_distances.jl
    julia -p 4 benchmark_distances.jl
    exit
}
