params.outdir = 'results'
process FILTER_OVERLAP_EXONS {
	tag "Filters out entries in which overlap with overlapping exons is detected, as well as overlap between exons due to divergent transcripts."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(bed), path(repeats)
    val exon_gap
    val cores
    val repeat_filter

    output:
    path "*.filtered.bed"
    path "*.LINE1.bed"

    script:

    """
    POLARS_MAX_THREADS=${cores} overlapping_exon_filter.py "${bed}" "${bam}" "${repeats}" -g ${exon_gap} -r ${repeat_filter}
    """
}