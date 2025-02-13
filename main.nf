#!/usr/bin/env nextflow

/*
* Pipeline for computing genome-wide coverage from fastq.gz files.
* Author:
* - Robert Linder <rlinder02@gmail.com>
*/

nextflow.enable.dsl = 2

/*
* Default pipeline parameters. They can be overriden on the command line eg.
* given `params.foo` specify on the run command line `--foo some_value`.
*/

params.species = 'human'
params.fastq = "$projectDir/full_assets/*.{fastq,fq}.gz"
params.cores = 20
params.align_cores = 8
params.minimap2_default = 'map-hifi'
params.bamtofastq = false
// set the maximum distance between spliced out exons
params.distance = 2500000
// set the fraction of an exon a read needs to overlap to continue on through the pipeline
params.overlap = 0.1
// set the fraction of a known repeat the exon-overlapped portion of a read needs to overlap to be filtered out
params.repeat = 0.2
params.repeat_spliced = 0.1
params.extra_gffs = false
// set the fraction of the intron length that needs to be spliced out for a read to continue through the pipeline
params.exon_gap = 0.2

/* 
* Load modules
*/

include {BAM_TO_FASTQ} from './modules/0_BAM_TO_FASTQ/main'
include {REFORMAT_ANNS} from './modules/1_REFORMAT_ANNS/main'
include {COMBINE_ANNS} from './modules/2_COMBINE_ANNS/main'
include {ALIGN} from './modules/3_ALIGN/main'
include {ALIGN_SPLICE_FULL} from './modules/4_ALIGN_SPLICE_FULL/main'
include {ONLY_SPLICE_DEL} from './modules/5_ONLY_SPLICE_DEL/main'
include {SNIFFLES} from './modules/6_SNIFFLES/main'
include {VCF_CONCAT} from './modules/7_VCF_CONCAT/main'
include {BAM_MERGE} from './modules/8_BAM_MERGE/main'
include {INS_FASTA} from './modules/9_INS_FASTA/main'
include {ALIGN_SPLICE} from './modules/10_ALIGN_SPLICE/main'
include {ONLY_SPLICED} from './modules/11_ONLY_SPLICED/main'
include {FILTER_DUPS} from './modules/12_FILTER_DUPS/main'
include {MARK_ALIGNMENTS} from './modules/13_MARK_ALIGNMENTS/main'
include {OVERLAP_FEATURES_INS} from './modules/14_OVERLAP_FEATURES_INS/main'
include {OVERLAP_LEN} from './modules/15_OVERLAP_LEN/main'
include {REFORMAT_BED} from './modules/16_REFORMAT_BED/main'
include {EXTRACT_READS} from './modules/17_EXTRACT_READS/main'
include {BAM_INTERSECT} from './modules/18_BAM_INTERSECT/main'
include {PS_FILTER} from './modules/19_PS_FILTER/main'
include {FILTER_BAM} from './modules/20_FILTER_BAM/main'
include {FILTER_OVERLAP_EXONS} from './modules/21_FILTER_OVERLAP_EXONS/main'
include {EXTRACT_READS_FINAL} from './modules/22_EXTRACT_READS_FINAL/main'
include {FIND_INS_READ} from './modules/23_FIND_INS_READ/main'
include {FILTER_FASTA} from './modules/24_FILTER_FASTA/main'
include {EXTRACT_FLANKS} from './modules/25_EXTRACT_FLANKS/main'
include {DELINEATE_INS} from './modules/26_DELINEATE_INS/main'
include {COMBINE_INS_REPORTS} from './modules/27_COMBINE_INS_REPORTS/main'

/* 
* main script workflow
*/

workflow {
    if (params.bamtofastq) {
        bamfiles = "$projectDir/full_assets/*.bam{,.pbi}"
        bams = Channel.fromFilePairs(bamfiles, checkIfExists: true)
        BAM_TO_FASTQ(bams, params.cores)
        fastq_ch = BAM_TO_FASTQ.out
        fastq_ch.view()
    } else {
        fastq_ch = Channel.fromPath(params.fastq)
    }

    if (params.species == 'human') {
            reference = "$projectDir/assets/GCF_009914755.1_T2T-CHM13v2.0_genomic_converted.fasta{,.fai}"
            primary_annotation = "$projectDir/assets/T2T-CHM13v2_genomic_converted_exons3.gff"
    } else if (params.species == 'mouse') {
            //reference = "$projectDir/assets/Mus_musculus.GRCm39.dna_sm.toplevel.corrected.fa{,.fai}"
            //primary_annotation = "$projectDir/assets/Mus_musculus.GRCm39.110_converted.gff"
            reference = "$projectDir/assets/BALB_cJ_v3_converted.fasta{,.fai}"
            primary_annotation = "$projectDir/assets/Mus_musculus-GCA_921997145.2-2023_06-genes_converted_ps_fixed.gff3"
    }
    
    if (params.extra_gffs) {
        pseudo_gffs = "$projectDir/full_assets/*.gff" 
        pseudo_gff_ch = Channel.fromPath(pseudo_gffs)
        sep_pseudo_gff_ch = pseudo_gff_ch.flatten()
        gff_ch1 = Channel.fromPath(primary_annotation)
        REFORMAT_ANNS(sep_pseudo_gff_ch, params.species)
        reformatted_gff_ch = REFORMAT_ANNS.out.collect()
        COMBINE_ANNS(gff_ch1, reformatted_gff_ch)
        gff_ch = COMBINE_ANNS.out
    } else {
        gff_ch = Channel.fromPath(primary_annotation)
    }
    ref_ch = Channel.fromFilePairs(reference)
    fastq_id = fastq_ch.map { file -> 
                def key = file.name.toString().split('/').last().split('.fastq.gz').first().split('.fq.gz').first()
                return tuple(key, file) }.groupTuple()
    combined_ch = fastq_id.combine(ref_ch)
    ALIGN(combined_ch, params.align_cores, params.minimap2_default)
    ALIGN_SPLICE_FULL(combined_ch, params.cores)
    ONLY_SPLICE_DEL(ALIGN_SPLICE_FULL.out[0], ALIGN_SPLICE_FULL.out[1])
    bam_id = ALIGN.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    bai_id = ALIGN.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    bam_splice_id = ONLY_SPLICE_DEL.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    bai_splice_id = ONLY_SPLICE_DEL.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    bam_bai_ch = bam_id.combine(bai_id, by:0)
    bam_bai_splice_ch = bam_splice_id.combine(bai_splice_id, by:0)
    ref_combined_full_ch = bam_bai_ch.combine(ref_ch)
    ref_combined_full_splice_ch = bam_bai_splice_ch.combine(ref_ch)
    all_ref_combined_full_ch = ref_combined_full_ch.concat(ref_combined_full_splice_ch)
    all_ref_combined_full_ch.view()
    SNIFFLES(all_ref_combined_full_ch, params.cores)
    vcf_id = SNIFFLES.out[0].collect().flatten().map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    vcf_idx_id = SNIFFLES.out[1].collect().flatten().map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    vcf_idx_ch = vcf_id.combine(vcf_idx_id, by:0)
    VCF_CONCAT(vcf_idx_ch)
    regular_bams_id = ALIGN.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_bais_id = ALIGN.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    spliced_bams_id = ALIGN_SPLICE_FULL.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    spliced_bais_id = ALIGN_SPLICE_FULL.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_bam_bai_combined_ch = regular_bams_id.combine(regular_bais_id, by:0)
    spliced_bam_bai_combined_ch = spliced_bams_id.combine(spliced_bais_id, by:0)
    all_bam_bai_combined_ch = regular_bam_bai_combined_ch.combine(spliced_bam_bai_combined_ch, by:0)
    BAM_MERGE(all_bam_bai_combined_ch, params.cores)
    INS_FASTA(VCF_CONCAT.out[0])
    fasta_id = INS_FASTA.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('.fasta').first()
                    return tuple(key, file) }.groupTuple()
    fasta_ref_ch = fasta_id.combine(ref_ch)
    ALIGN_SPLICE(fasta_ref_ch, params.cores, params.distance)
    ONLY_SPLICED(ALIGN_SPLICE.out[0], ALIGN_SPLICE.out[1])
    FILTER_DUPS(ONLY_SPLICED.out[0], ONLY_SPLICED.out[1])
    filtered_dups_bam_ch = FILTER_DUPS.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    filtered_dups_bai_ch = FILTER_DUPS.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    filtered_dups_ch = filtered_dups_bam_ch.combine(filtered_dups_bai_ch, by:0)
    MARK_ALIGNMENTS(filtered_dups_ch)
    regular_filtered_bam_ch = MARK_ALIGNMENTS.out.filter( ~/.*insertions.*/ )
    marked_bam_id = MARK_ALIGNMENTS.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('.bam').first()
                    return tuple(key, file) }.groupTuple()
    ann_combined_bam_ch = marked_bam_id.combine(gff_ch)
    OVERLAP_FEATURES_INS(ann_combined_bam_ch)
    OVERLAP_LEN(OVERLAP_FEATURES_INS.out[0], params.overlap)
    ins_combined_ch = OVERLAP_LEN.out.combine(gff_ch) 
    REFORMAT_BED(ins_combined_ch)
    regular_reformatted_bed_ch = REFORMAT_BED.out.filter( ~/.*insertions.*/ )
    vcf_id2 = VCF_CONCAT.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    genc_id = regular_reformatted_bed_ch.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    vcf_genc_id_ch = vcf_id2.combine(genc_id, by:0)
    merged_bam_id = BAM_MERGE.out[0].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    merged_bai_id = BAM_MERGE.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    merged_bam_bai_ch = merged_bam_id.combine(merged_bai_id, by:0)
    merged_bam_bai_vcf_genc_id_ch = merged_bam_bai_ch.combine(vcf_genc_id_ch, by:0)
    EXTRACT_READS(merged_bam_bai_vcf_genc_id_ch)
    extract_bam_id = EXTRACT_READS.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    bam_bed_ch = genc_id.combine(extract_bam_id, by:0)
    ann_bam_bed_ch = bam_bed_ch.combine(gff_ch)
    BAM_INTERSECT(ann_bam_bed_ch)
    PS_FILTER(BAM_INTERSECT.out, params.cores)
    regular_ps_filter_ch = PS_FILTER.out[0].filter( ~/.*annotated.all.*/)
    ps_regular_filter_id = regular_ps_filter_ch.map { file ->
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    ps_regular_bam_ch = extract_bam_id.combine(ps_regular_filter_id, by:0)
    FILTER_BAM(ps_regular_bam_ch)
    regular_pseudogneless_bed_ch = PS_FILTER.out[1].filter( ~/.*annotated.all.*/ )

    // Need a channel for the bed file with overlapping repeats 
    regular_repeat_genc_id = OVERLAP_FEATURES_INS.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple() 
    regular_ps_genc_id = regular_pseudogneless_bed_ch.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_ps_bam_id = regular_filtered_bam_ch.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_ps_genc_bam_ch = regular_ps_bam_id.combine(regular_ps_genc_id, by:0)
    regular_ps_repeat_genc_bam_ch = regular_ps_genc_bam_ch.combine(regular_repeat_genc_id, by:0)
    FILTER_OVERLAP_EXONS(regular_ps_repeat_genc_bam_ch, params.exon_gap, params.cores, params.repeat)
    regular_filtered_bed_ch = FILTER_OVERLAP_EXONS.out[0].filter( ~/.*annotated.all.*/ )    
    regular_filtered_bam_ch2 = FILTER_BAM.out[0].filter( ~/.*merged.*/ )
    regular_genc_id = regular_filtered_bed_ch.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_bam_id = regular_filtered_bam_ch2.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_genc_bam_ch = regular_bam_id.combine(regular_genc_id, by:0)
    EXTRACT_READS_FINAL(regular_genc_bam_ch)
    regular_filtered_bam_ch3 = EXTRACT_READS_FINAL.out.filter( ~/.*merged.*/ )
    ins_location_id = VCF_CONCAT.out[2].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_filtered_bam_id3 = regular_filtered_bam_ch3.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    ins_location_filtered_bam_id3_ch = ins_location_id.combine(regular_filtered_bam_id3, by:0)
    FIND_INS_READ(ins_location_filtered_bam_id3_ch)
    fasta_id_new = INS_FASTA.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_ps_id_ch = PS_FILTER.out[2].filter( ~/.*annotated.all.*/)
    no_ps_filter_id = regular_ps_id_ch.map { file ->
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    no_ps_fasta_ch = fasta_id_new.combine(no_ps_filter_id, by:0)
    FILTER_FASTA(no_ps_fasta_ch)
    filtered_fasta_id = FILTER_FASTA.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    filtered_bam_id = FIND_INS_READ.out.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    ins_loc_fasta_ch = ins_location_id.combine(filtered_fasta_id, by:0)
    ins_loc_filtered_bam_ch = ins_loc_fasta_ch.combine(filtered_bam_id, by:0)
    EXTRACT_FLANKS(ins_loc_filtered_bam_ch)
    reg_filtered_bam_id = regular_filtered_bam_ch.map { file ->
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    regular_genc_bed_bam_ch = regular_genc_id.combine(reg_filtered_bam_id, by:0)
    DELINEATE_INS(regular_genc_bed_bam_ch)
    ins_report_regular = DELINEATE_INS.out[0].filter( ~/.*ins.*/ )
    ins_report_ch = ins_report_regular.map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    flank_report_ch = EXTRACT_FLANKS.out[1].map { file -> 
                    def key = file.name.toString().split('/').last().split('\\.').first()
                    return tuple(key, file) }.groupTuple()
    ins_flank_ch = ins_report_ch.combine(flank_report_ch, by:0)
    COMBINE_INS_REPORTS(ins_flank_ch)
} 