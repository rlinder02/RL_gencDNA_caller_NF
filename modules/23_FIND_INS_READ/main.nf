params.outdir = 'results'
process FIND_INS_READ {
	tag "If there are multiple alignments from the two separate mapping steps that map to the same insertion, find the alignment in which the insertion was most likely found."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(ins_locs), path(bam)

    output:
    path "*.ins.reads.bam"

    script:

    """
    extract_ins_reads.py "${ins_locs}" "${bam}"
    """
}