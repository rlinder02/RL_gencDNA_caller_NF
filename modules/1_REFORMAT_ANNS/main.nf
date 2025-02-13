params.outdir = 'results'
process REFORMAT_ANNS {
	tag "Reformat additional annotation gff files so they are compatible with downstream processes" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path gff
    val species 

    output:
    path "*.additional.gff"
    
    script:

    """
    combine_gff_files.py "${gff}" "${species}"
    """
}