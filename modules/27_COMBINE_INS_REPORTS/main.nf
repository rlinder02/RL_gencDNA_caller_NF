params.outdir = 'results'
process COMBINE_INS_REPORTS {
	tag "Combine the flanking and insertion reports into a single file, merging by insertion id" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(ins), path(flanks) 

    output:
    path "*.ins_flank.report.txt"

    script:

    """
    merge_ins_flank_reports.py "${ins}" "${flanks}" -o "${sample_id}"
    """
}