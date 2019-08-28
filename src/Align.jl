function minimap2consensus(
    ref_fasta_path, query, out_bam_path, out_fasta_path;
    preset, min_ac, min_af, max_dp, threads=3)

    run(pipeline(`minimap2 -ax $preset -t $threads $ref_fasta_path $query`,
                 `samtools sort -O bam -o $out_bam_path -`))

    run(`samtools index $out_bam_path`)

    open(out_fasta_path, "w") do out_fasta_stream
        bam2pileup2consensus(out_bam_path, out_fasta_stream, ref_fasta_path,
                             min_ac, min_af, max_dp)
    end
end

"""
    align_short_reads(ref_fasta_path, fastq_list, out_bam_path, out_fasta_path; <keyword arguments>)

Align short reads against a reference, and produce a consensus FASTA.

Reads are aligned using `minimap2` with preset "sr". To generate a
consensus FASTA, a pileup is generated with `samtools mpileup`, and
the consensus allele is output at each site passing filters (see
keyword arguments); the missing allele "-" is output at sites not
passing filters. Note the consensus FASTA does not contain insertions,
only substitutions and deletions.

# Arguments
- `ref_fasta_path`: filename of reference FASTA. Should be gzipped and
   indexed with `samtools faidx`.
- `fastq_list`: list of fastq files for the sample.
- `out_bam_path`: Output BAM filename.
- `out_fasta_path`: Output consensus FASTA filename.
- `min_ac=10`: Filter for minimum count of the consensus allele.
- `min_af=0.9`: Filter for minimum frequency of the consensus allele.
- `max_dp=typemax(Int)`: Filter for maximum read depth to consider.
- `threads=3`: Number of threads for minimap2 (Note: pileup stage not currently parallelized)
"""
function align_short_reads(
    ref_fasta_path, fastq_list, out_bam_path, out_fasta_path;
    min_ac=10, min_af=0.9, max_dp=typemax(Int), threads=3)

    minimap2consensus(
        ref_fasta_path, fastq_list, out_bam_path, out_fasta_path,
        preset="sr", min_ac=min_ac, min_af=min_af, max_dp=max_dp,
        threads=threads)
end

"""
    align_assembly(ref_fasta_path, asm_fasta_path, out_bam_path, out_fasta_path; <keyword arguments>)

Align assembled FASTA against a reference.

Reads are aligned using `minimap2` (see keyword arguments for
settings). Then, a pileup is generated with `samtools mpileup`, and
the consensus allele is output at each site passing filters (see
keyword arguments); the missing allele "-" is output at sites not
passing filters. Note the consensus FASTA does not contain insertions,
only substitutions and deletions. Also note that the consensus FASTA
uses the contig names from the reference FASTA, not the new assembly.

# Arguments
- `ref_fasta_path`: filename of reference FASTA. Should be gzipped and
   indexed with `samtools faidx`.
- `asm_fasta_path`: filename of assembly FASTA.
- `out_bam_path`: Output BAM filename.
- `out_fasta_path`: Output consensus FASTA filename.
- `preset="asm5"`: Preset for minimap2. asm5 means <=5% divergence. Other relevant presets are asm10 and asm20.
- `min_ac=1`: Filter for minimum count of the consensus allele.
- `min_af=1`: Filter for minimum frequency of the consensus allele.
- `max_dp=1`: Filter for maximum read depth to consider.
- `threads=3`: Number of threads for minimap2 (Note: pileup stage not currently parallelized)
"""
function align_assembly(
    ref_fasta_path, asm_fasta_path, out_bam_path, out_fasta_path;
    preset="asm5", min_ac=1, min_af=1, max_dp=1, threads=3)

    minimap2consensus(
        ref_fasta_path, asm_fasta_path, out_bam_path, out_fasta_path,
        preset=preset, min_ac=min_ac, min_af=min_af, max_dp=max_dp,
        threads=threads)
end

export align_short_reads
export align_assembly
