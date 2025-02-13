# RL_gencDNA_caller_NF

This is a Nextflow pipeline for analyzing WGS data generated from Pacific Biosciences HiFi reads. The specific aim is to find evidence for the presence of cDNA that has retro-inserted into the genome via reverse trasncription of processed mRNA coupled with double-stranded DNA damage. These genomic features are known as gencDNA (genomic cDNA).

[![homepage](https://img.shields.io/badge/nextflow-%E2%89%A523.10.0-brightgreen.svg)](https://nextflow.io/ "Redirect to nextflow homepage")

## Requirements

- Unix-like operating system (Linux, macOS, etc)
- Java 8 or later
- Docker 24.0.6 or later
- Nextflow 23.10.0 or later

## Quickstart

If you don't already have Nextflow, install by using the following command:

```
curl -s https://get.nextflow.io | bash
```

The Dockerfile and conda.yaml are included if updated or additional software tools are desired to be added to the Docker image used throughout this pipeline. To download the Docker image from dockerhub, use this command:

```
docker pull rlinder02/gencdna_indel:v1.0.11
```

Clone the pipeline and then run on your local machine from GitHub:

```
git clone https://github.com/rlinder02/RL_gencDNA_caller_NF.git
```

To run this pipeline, you will need PacBio HiFi reads in either `.fastq.gz` or `.bam` format (with the accompanying `.pbi` file for bam files) in the `full_assets` directory. If these are bam files, they will first be converted to `.fastq.gz` files before proceeding with the pipeline as Minimap2 does not currently accept bam files as input by adding `--bamtofastq` to the command line argument. For analyzing BALB mouse samples, you will need to include the `BALB_cJ_v3_converted.fasta`, `BALB_cJ_v3_converted.fasta.fai`, and `Mus_musculus-GCA_921997145.2-2023_06-genes_converted_ps_fixed.gff3` files in the `assets` folder. The `BALB_cJ_v3_repeats.gff`, `Mus_musculus-GCA_921997145.2-2023_06-genes_converted_ps_fixed.gff3`, and `Mus_musculus.GRCm39_to_Balb_repeats_converted.fixed.gff` files must be included in the `full_assets` folder. These files can be found in my S3 bucket in the `gencDNA/gencDNA_caller_annotation_files/mouse` directory. These are the default files for when the species is set as mouse.

Alternatively, for analyzing other mouse strains, include the `Mus_musculus.GRCm39.dna_sm.toplevel.corrected.fa`, `Mus_musculus.GRCm39.dna_sm.toplevel.corrected.fa.fai`, and `Mus_musculus.GRCm39.110_converted.gff` files in the `assets` folder. The `Mus_musculus.GRCm39_repeats_converted.gff` file must be included in the `full_assets` folder. These files can be found in my S3 bucket in the `gencDNA/gencDNA_caller_annotation_files/mouse_C57` directory. If using these files, go into the `main.nf` file and comment out the default mouse reference and primary annotation variables, while removing the comments from the alternative variables just above. 

Analyzing mouse samples involves adding `--species mouse` to the command line argument. 

For analyzing human data, include the `GCF_009914755.1_T2T-CHM13v2.0_genomic_converted.fasta`, `GCF_009914755.1_T2T-CHM13v2.0_genomic_converted.fasta.fai`, and `T2T-CHM13v2_genomic_converted_exons3.gff` files in the `assets` folder. The `GCF_009914755.1_T2T-CHM13v2.0_repeats.gff`, `gencode.v45.hg38_to_T2T_ps_annotation.gff`, and `GRCh38.p14_to_T2T_repeats_converted.fixed.gff` files must be included in the `full_assets` folder. These files can be found in my S3 bucket in the `gencDNA/gencDNA_caller_annotation_files/human` directory.


Once these files are in the `assets` and `full_assets` directory, run the following code snippet (assumes the species is human and the input files are in `.fastq.gz` format), replacing the default number of cores and align_cores with the numbers you want. Be careful if running multiple samples in parallel, as the total number of CPUs used will simultaneously will be the sum of these number multiplied by the number of samples. The default output is a text file with information on the potential gencDNAs found:

```
nextflow run main.nf --extra_gffs  
```
To run in the background, run:

```
nextflow run -bg main.nf --extra_gffs
```

After the pipeline has completed, the text files should be in the `results` folder. At this point, if you no longer need the intermediate files, delete the contents of the `full_assets` folder and the `work` directory produced during the run to free up space on the drive.

## Pipeline parameters

At the command line, you can specify several options (documented [here](https://www.nextflow.io/docs/latest/)). The most relevant options for this pipeline are listed as `params.` Variables are towards the top the main.nf file and can be specified on the command line line when launching Nextflow.

## Pipeline results

The results are copied to a generic `results` folder.

## Credits

The general format of the pipeline was inspired by the work done [here](https://github.com/nextflow-io/rnaseq-nf) and [here](https://github.com/nf-core/rnaseq).
