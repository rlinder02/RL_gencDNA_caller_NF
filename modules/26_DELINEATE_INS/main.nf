params.outdir = 'results'
process DELINEATE_INS {
	tag "Extract insertion sequences from the final list of insertions that passed all previous filters" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bed), path(bam) 

    output:
    path "*.report.txt"
    path "*.bai"

    script:

    """
    samtools index "${bam}"
    delineate_aligned_insertions.py "${bam}" "${bed}"
    """
}