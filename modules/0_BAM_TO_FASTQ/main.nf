params.outdir = 'results'

process BAM_TO_FASTQ {
	tag "Convert bam files to fastq.gz files; needs an index file to enable conversion"
	label 'preprocess'
	publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam)
    val cores  

    output:
    path "*.fastq.gz" 

    shell:
    '''
    bam2fastq -o !{sample_id} !{bam[0]} -j !{cores}
    '''
}