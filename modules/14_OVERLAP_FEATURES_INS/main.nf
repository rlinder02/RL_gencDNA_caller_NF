params.outdir = 'results'
process OVERLAP_FEATURES_INS {
	tag "Find overlapping features in insertions; looking for overlapped exons of any genes- can specify the fraction of the feature that must be covered to call" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(bam), path(all_ann)

    output:
    path "*.annotated.exons.bed"
    path "*.repeats.annotated.bed"

    shell:

    '''
    cat !{all_ann} | awk '$2!="RepeatMasker"' > no_repeats.gff
    bedtools intersect -split -wo -bed -abam !{bam} -b no_repeats.gff > !{sample_id}.annotated.bed
    cat !{sample_id}.annotated.bed | awk '$15=="exon"' > !{sample_id}.annotated.exons.bed 
    cat !{all_ann} | awk '$2=="RepeatMasker"' > only_repeats.gff
    bedtools intersect -split -wo -bed -abam !{bam} -b only_repeats.gff > !{sample_id}.repeats.annotated.bed   
    '''
}