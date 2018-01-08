# pipeline_rotation
A repository containing my pipeline which was developed as part of the CBI1 rotation (Nov - Dec 2017).
Briefly, the pipeline consists of the following stages:
* QC on raw, unprocessed reads using FastQC (v0.11.5)
* Aligning short reads to hg19 reference genome using BWA-MEM (v0.7.15)
* QC step using Samtools (v1.3.1)
* Indel realignment using GATK (v3.6)
* Mark duplicates using Sambamba (v0.6.3)
* Read group information added using Picard Tools (v2.5.0)
* Coverage was calculated using Sambamba (v0.6.3)
* On/off target reads were calculated using Picard Tools (v2.5.0)
* Variant calling using GATK (v3.6)
* Variant annotation using ExAC (common_all_20161122.vcf.gz) and SnpSift (v4.3p)



# Usage
The pipeline can be run using the following commands in a bash terminal:

`python pipeline.py worklist-sampleID broad_bed small_bed`

A read world example would be:

`python pipeline.py 1606034-S1612259 /example_fastqs/beds/CTDFinaldesignwith25bp_v3.bed /example_fastqs/beds/CTD_EDS_C_25_v4.bed`

