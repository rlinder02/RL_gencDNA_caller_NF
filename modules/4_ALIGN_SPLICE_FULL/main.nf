params.outdir = 'results'
process ALIGN_SPLICE_FULL {
	tag "Align all reads to a reference genome in splice-aware mode and pull out alignments with large deletions and/or spliced-out introns (Ns), as well as running the resulting bam file through Sniffles2 again to detect insertions that occur adjacent to deletions"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(fastq), val(ref_id), path(ref)
    val cores   

    output:
    path "*.sorted.bam" 
    path "*.sorted.bam.bai"

    shell:

    '''
    zcat !{fastq} | sed -E 's/(@m[0-9][0-9][0-9][0-9][0-9]_)/\\1SPLICE_/' > !{sample_id}.splice.fastq
    minimap2 --eqx -t !{cores} -Yax splice:hq -G 500000 -o !{sample_id}.spliceaware.sam !{ref[0]} !{sample_id}.splice.fastq
    rm !{sample_id}.splice.fastq
    samtools sort -m 4G -o !{sample_id}.spliceaware.sorted.bam -@ !{cores} !{sample_id}.spliceaware.sam
    samtools index -@ !{cores} !{sample_id}.spliceaware.sorted.bam
    rm !{sample_id}.spliceaware.sam
    '''
}