params.outdir = 'results'
process PS_FILTER {
	tag "Extracts reads that overlap known pseudogenes and purges the bed file of these reads."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path bed
    val cores

    output:
    path "*.pseudogene.reads.txt"
    path "*.pseudogeneless.bed"
    path "*.pseudogeneless.reads.txt"

    script:

    """
    POLARS_MAX_THREADS=${cores} pseudogene_readnames.py "${bed}"
    """
}