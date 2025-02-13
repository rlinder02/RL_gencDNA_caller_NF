params.outdir = 'results'
process REFORMAT_BED {
	tag "Reformats the bed files to add info for further processing."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple path(bed), path(gff)

    output:
    path "*.reformatted.bed"

    script:

    """
    reformat_bed.py "${gff}" "${bed}"
    """
}