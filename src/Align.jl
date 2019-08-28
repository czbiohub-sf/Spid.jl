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

function align_short_reads(
    ref_fasta_path, fastq_list, out_bam_path, out_fasta_path;
    min_ac=10, min_af=0.9, max_dp=typemax(Int), threads=3)

    minimap2consensus(
        ref_fasta_path, fastq_list, out_bam_path, out_fasta_path,
        preset="sr", min_ac=min_ac, min_af=min_af, max_dp=max_dp,
        threads=threads)
end

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
