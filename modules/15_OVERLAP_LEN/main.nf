params.outdir = 'results'
process OVERLAP_LEN {
	tag "Filters by the amount of overlap between the insertion and overlapping exon as well as the insertion and any overlapping repetitive elements."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path bed
    val overlap

    output:
    path "*.overlap.bed"

    script:

    """
    exon_overlap_filter.py "${bed}" -o "${overlap}"
    """
}