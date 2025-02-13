params.outdir = 'results'
process VCF_CONCAT {
	tag "Take the vcf files generated from the normal and splice-aware Sniffles2 insertion calls and combine into one."
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(vcfs), path(idxs)

    output:
    path "*.seqs.txt" 
    path "*insertions.vcf"
    path "*.locations.txt"

    shell:

    '''
    bcftools concat -a -D -o !{sample_id}.vcf !{vcfs[0]} !{vcfs[1]}
    cat !{sample_id}.vcf | grep "Sniffles2.INS" | cut -f1,2,3,5,6,8 | awk -F'\t' -v OFS='\t' '{$3=$3 "." (++count[$3])}1' > !{sample_id}.insertions.vcf
    cat !{sample_id}.insertions.vcf | cut -f3,4 > !{sample_id}.insertions.seqs.txt
    cat !{sample_id}.insertions.vcf | cut -f1,2,3,6 > !{sample_id}.insertions.locations.txt
    '''
}