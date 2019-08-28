FROM debian:unstable-20190812

ENV JULIA_DEPOT_PATH=/opt/julia
RUN apt-get update && \
        apt-get install -y julia samtools minimap2 curl && \
        mkdir -p $JULIA_DEPOT_PATH

ENV JULIA_PROJECT=/opt/Spid.jl
COPY src/*.jl $JULIA_PROJECT/src/
COPY Project.toml $JULIA_PROJECT/
COPY bin/spid.jl /usr/local/bin/
RUN julia -O3 -e 'using Pkg; Pkg.instantiate(); using Spid' && \
        # change permissions AFTER installing julia packages
        chmod -R a+rw /opt/julia
