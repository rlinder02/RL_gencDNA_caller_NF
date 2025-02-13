params.outdir = 'results'
process SNIFFLES {
	tag "Take the sorted and indexed bam files and use them to call SVs using Sniffles2."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(bai), val(ref_id), path(ref)
    val cores  

    output:
    path "*.combined.vcf.gz"
    path "*.combined.vcf.gz.tbi"

    shell:

    '''
    sniffles \\
    --input !{bam} \\
    --vcf !{sample_id}.vcf.gz \\
    --reference !{ref[0]} \\
    --output-rnames \\
    --mosaic \\
    --mosaic-af-min 0 \\
    --minsupport 1 \\
    --mapq 0 \\
    --threads !{cores}
    sniffles \\
    --input !{bam} \\
    --vcf !{sample_id}2.vcf.gz \\
    --reference !{ref[0]} \\
    --output-rnames \\
    --minsupport 1 \\
    --mapq 0 \\
    --threads !{cores}
    bcftools concat -a -D -o !{sample_id}.combined.vcf.gz !{sample_id}.vcf.gz !{sample_id}2.vcf.gz
    bcftools index -t !{sample_id}.combined.vcf.gz
    '''
}