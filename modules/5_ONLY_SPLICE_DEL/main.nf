params.outdir = 'results'
process ONLY_SPLICE_DEL {
	tag "Filter out reads that do not have an N or a 30+D in the cigar string"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path bam
    path bai

    output:
    path "*deletions.bam"
    path "*deletions.bam.bai" 

    script:

    """
    output_spliced_alignments_only.py "${bam}"
    """
}