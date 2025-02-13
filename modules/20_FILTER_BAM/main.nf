params.outdir = 'results'
process FILTER_BAM {
	tag "Filter out reads that overlap known pseudogenes from the original aligned bam file."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(ps_reads)

    output:
    path "*.bam"

    script:

    """
    filter_out_pseudogenes.py "${bam}" "${ps_reads}"
    """
}