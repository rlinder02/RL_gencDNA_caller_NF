params.outdir = 'results'
process ALIGN {
	tag "Align fastq.gz files (with read names modified to reflect mapping type) to a reference genome (only human and mouse currently supported), sort, and index"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(fastq), val(ref_id), path(ref)
    val cores
    val mm2_default

    output:
    path "*.bam" 
    path "*.bai"

    shell:

    '''
    zcat !{fastq} | sed -E 's/(@m[0-9][0-9][0-9][0-9][0-9]_)/\\1REG_/' > !{sample_id}.reg.fastq
    minimap2 --eqx --secondary=no -t !{cores} -Yax !{mm2_default} -o !{sample_id}.sam !{ref[0]} !{sample_id}.reg.fastq
    rm !{sample_id}.reg.fastq
    samtools sort -m 4G -o !{sample_id}.sorted.bam -@ !{cores} !{sample_id}.sam
    samtools index -@ !{cores} !{sample_id}.sorted.bam
    rm !{sample_id}.sam
    '''
}