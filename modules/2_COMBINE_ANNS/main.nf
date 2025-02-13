params.outdir = 'results'
process COMBINE_ANNS {
	tag "Combine the reformatted gffs with the reference gff" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
     path ref_gff
     path ps_gffs

    output:
    path "*.combinedannotations.gff"

    shell:

    '''
    gff_name="$(basename  !{ref_gff} .gff)".combinedannotations.gff
    cat !{ref_gff} !{ps_gffs} > $gff_name
    '''
}