params.outdir = 'results'
process INS_FASTA {
	tag "Convert the file of insertion ids and sequences to fasta format for alignment"
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    path(ins_file)

    output:
    path "*.fasta" 
    
    script:
    
    """
    mapping_insertions.py "${ins_file}"
    """
}