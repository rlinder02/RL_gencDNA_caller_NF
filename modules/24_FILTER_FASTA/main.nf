params.outdir = 'results'
process FILTER_FASTA {
	tag "Extract insertion sequences from the final list of insertions that passed all previous filters" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(fasta), path(ins_ids)

    output:
    path "*.filtered.fasta"

    script:

    """
    filter_insertion_fasta.py "${ins_ids}" "${fasta}"
    """
}