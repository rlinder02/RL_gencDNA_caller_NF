params.outdir = 'results'
process EXTRACT_READS {
	tag "Extracts reads with potential gencDNA insertions."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(bed)

    output:
    path "*.filtered.bam"

    script:

    """
    extract_reads_bamfile.py "${bed}" "${vcf}" "${bam}"
    """
}