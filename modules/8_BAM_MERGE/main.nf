params.outdir = 'results'
process BAM_MERGE {
	tag "Take the bam files generated from the normal and splice-aware Minimap2 run and combine into one."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(reg_bam), path(reg_bai), path(splice_bam), path(splice_bai)
    val cores

    output:
    path "*.merged.bam" 
    path "*.merged.bam.bai"

    shell:

    '''
    samtools merge -c -@ !{cores} -o !{sample_id}.merged.bam !{reg_bam} !{splice_bam}
    samtools index !{sample_id}.merged.bam
    '''
}