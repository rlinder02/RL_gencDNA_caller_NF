params.outdir = 'results'
process ALIGN_SPLICE {
	tag "Align insertion fasta sequences to a reference genome in splice-aware mode and pull out alignments with large deletions and/or spliced-out introns (Ns)"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(fasta), val(ref_id), path(ref)
    val cores
    val distance   

    output:
    path "*.sorted.bam" 
    path "*.sorted.bam.bai"

    shell:

    '''
    minimap2 --eqx -t !{cores} -Yax splice:hq -G !{distance} -o !{sample_id}.sam !{ref[0]} !{fasta}
    samtools sort -m 4G -o !{sample_id}.sorted.bam -@ !{cores} !{sample_id}.sam
    samtools index -@ !{cores} !{sample_id}.sorted.bam
    '''
}