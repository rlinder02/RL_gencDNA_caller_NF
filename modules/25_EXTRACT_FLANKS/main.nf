params.outdir = 'results'
process EXTRACT_FLANKS {
	tag "Extract sequences flanking insertions" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(ins_ids), path(fasta), path(bam)

    output:
    path "*.fasta"
    path "*.report.txt"

    script:

    """
    extract_flanks.py "${fasta}" "${ins_ids}" "${bam}"
    """
}