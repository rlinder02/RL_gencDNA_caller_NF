params.outdir = 'results'
process BAM_INTERSECT {
	tag "Find overlapping features in the original bam (this is to filter out hits on known pseudogenes); don't split by introns, as this will at times miss overlaps with pseudogenes" 
	label 'calling'
    publishDir params.outdir, mode:'copy'
	debug true

    input:
    tuple val(sample_id), path(genc), path(bam), path(all_ann)

    output:
    path "*.annotated.all.fixed.tab.bed"

    shell:

    '''
    bedtools intersect -wo -bed -abam !{bam} -b !{all_ann} > !{sample_id}.annotated.bed
    cat !{sample_id}.annotated.bed | awk '$15!="region"' > !{sample_id}.annotated.noregion.bed
    cat !{sample_id}.annotated.noregion.bed !{genc} > !{sample_id}.annotated.all.tsv
    xsv fixlengths !{sample_id}.annotated.all.tsv > !{sample_id}.annotated.all.fixed.bed
    cat !{sample_id}.annotated.all.fixed.bed | awk -F'"' -v OFS='"' '{gsub(",", ".", $2);gsub(",", ".", $4);gsub(",", ".", $6);gsub(",", ".", $8)}1' | tr -d '"' | sed 's/,/\\t/g' > !{sample_id}.annotated.all.fixed.tab.bed
    '''
}