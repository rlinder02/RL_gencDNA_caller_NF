params.outdir = 'results'
process FILTER_DUPS {
	tag "Filter duplicate entries"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path bam
    path bai

    output:
    path "*nodups.bam"
    path "*nodups.bam.bai" 

    script:

    """
    filter_duplicate_entries.py "${bam}"
    """
}