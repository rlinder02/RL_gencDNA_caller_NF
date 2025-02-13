params.outdir = 'results'
process MARK_ALIGNMENTS {
	tag "Give a distinguishing mark to primary vs secondary alignments"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.marked.sorted.bam"

    script:

    """
    distinguishing_alignments.py "${bam}"
    """
}