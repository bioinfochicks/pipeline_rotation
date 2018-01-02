#activate virtual envs
#source /venvs/stp/bin/activate

import subprocess
import os
import getopt
import sys
import glob

#define some global variables
fastqc_path = "/software/FastQC/fastqc"
fastqc_out = "/home/stpuser/FastQC_output/"
sample_id = sys.argv[1]
broad_bed = sys.argv[2]
small_bed = sys.argv[3]
sample_list = glob.glob("/example_fastqs/" + sample_id + "*.fastq.gz")
out_dir = "/home/stpuser/Results/{id}".format(id = sample_id)
reference = "/reference_files/ucsc.hg19.nohap.masked.fasta"
dbsnp = "/reference_files/dbsnp/common_all_20161122.vcf.gz"

def run_cmd(cmd):
	try:
		subprocess.call(cmd, shell=True)
	except subprocess.CallProcessError as e: 
		log(e)
		print(e.output)

def make_dir(dir_path):
	os.mkdir(dir_path)
	print()
	return
	
def log(to_log):
	fname = "/home/stpuser/Results/{id}/log.txt".format(id = sample_id)
	with open(fname, "wb") as log_file:
		log_file.write(to_log)
		log_file.close()
		print("...logged.")
		
def run_fastqc(sample_list):
	#run fastqc on listed samples
	for file in sample_list:
		cmd = fastqc_path + " " + file + " " + "--outdir=" + fastqc_out 
		print(cmd)
		run_cmd(cmd)
	return 
			
def concat_files(sample_id, sample_list, out_dir):
	"""
	Gets all fastqs for sample and separates into each lane and read. Then, cats the L001 and L002 R1 files into R1.fastq.gz and L001 and L002 R2 files into R2.fastq.gz
	"""
	if os.path.isdir(out_dir):
		pass
	else:
		make_dir(out_dir)
		print("Made new directory for results")
		
	for file in sample_list:
		if "L001_R1" in file:
			L001_R1 = file
		elif "L002_R1" in file:
			L002_R1 = file
		elif "L001_R2" in file:
			L001_R2 = file
		elif "L002_R2" in file:
			L002_R2 = file

	cmd1 = "cat {0} {1} > {2}/R1.fastq.gz".format(L001_R1,L002_R1,out_dir)
	cmd2 = "cat {0} {1} > {2}/R2.fastq.gz".format(L001_R2,L002_R2,out_dir)
	
	run_cmd(cmd1)
	run_cmd(cmd2)
	return
	
def run_alignment(sample_id, sample_list):
	temp1 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_temp1.bam"
	temp2 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_temp2.bam"
	
	outdir = "/home/stpuser/Results/{id}".format(id = sample_id)
	temp3 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged.bam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted.bam"
	
	if os.path.isdir(out_dir):
		pass
	else:
		make_dir(out_dir)
		print("Made new directory for results")
		
	for file in sample_list:
		if "L001_R1" in file:
			L001_R1 = file
		elif "L002_R1" in file:
			L002_R1 = file
		elif "L001_R2" in file:
			L001_R2 = file
		elif "L002_R2" in file:
			L002_R2 = file
	
	cmd1 = "bwa mem " + reference + " " + L001_R1 + " " + L001_R2 + " | samtools view -bhS - > " + temp1 
	cmd2 = "bwa mem " + reference + " " + L002_R1 + " " + L002_R2 + " | samtools view -bhS - > " + temp2
	cmd3 = "samtools merge " + temp3 + " " + temp1 + " " + temp2
	cmd4 = "samtools sort " + temp3 + " -o " + out_bam + " | samtools index - "
	run_cmd(cmd1)
	run_cmd(cmd2)
	run_cmd(cmd3)

	os.remove(temp1)
	os.remove(temp2)
	os.remove(temp3)
	
	#https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups	
	
def index_bam(bam):
	cmd = "samtools index " + bam
	run_cmd(cmd)

def qc_read_mappings(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted.bam"
	out_text = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_samflagstats.txt"
	cmd = "samtools flagstat " + in_bam + " > " + out_text
	print("QC mapping cmd: " + cmd)
	run_cmd(cmd) 
	
def mark_dups(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted.bam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked.bam"
	cmd = "sambamba markdup -p " + in_bam + " " + out_bam
	print("Mark dups cmd: " + cmd)
	run_cmd(cmd)
	
def find_coverage(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG_indelsrealigned.bam"
	out1 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_coverage_depth_base_small_panel.txt"
	out2 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_coverage_depth_region_small_panel.txt"
	sorted_bed = "/home/stpuser/Results/" + sample_id + "/smallpanel_sorted.bed"
		
	cmd = "/software/bedtools2/bin/bedtools sort -i " + small_bed + " > " + sorted_bed
	run_cmd(cmd)
	cmd1 = "sambamba depth base -L " + sorted_bed + " --min-coverage=0 -q29 -m "  + in_bam + " > " + out1 
	print(cmd1)
	run_cmd(cmd1) 
	cmd2 = "sambamba depth region -L " + sorted_bed + " --cov-threshold=0 --cov-threshold=10 --cov-threshold=20 --cov-threshold=30 --cov-threshold=40 --cov-threshold=50 -q29 -m " + in_bam + " > " + out2
	print(cmd2)
	run_cmd(cmd2)

def add_readgroup(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked.bam"
	temp_bam = "/home/stpuser/Results/" + sample_id + "/temp.bam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG.bam"
	
	cmd = "java -jar /software/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=" + in_bam + " O=" + temp_bam + " RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=" + sample_id
	cmd2 = "samtools sort " + temp_bam + " -o " + out_bam + " | samtools index - "
	run_cmd(cmd)
	run_cmd(cmd2)
	
# def off_target_reads(sample_id):
def realign_indels(sample_id): 
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG.bam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG_indelsrealigned.bam"
	intervals = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_forIndelRealigner.intervals"
	regions = broad_bed

	index_bam(in_bam) # remove this after next run
	cmd1 = "java -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + reference + " -I " + in_bam + " -dt NONE -o " + intervals + " -L " + regions 
	cmd2 = "java -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R " + reference + " -I " + in_bam + " -targetIntervals " + intervals + " -dt NONE -o " + out_bam + " -L " + regions 
	run_cmd(cmd1) 
	run_cmd(cmd2)
	print(cmd1)
	print(cmd2)
	index_bam(out_bam)
	
def get_hs_stats(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG_indelsrealigned.bam"
	bait = "/home/stpuser/Results/" + sample_id + "/bait.intervals"
	intervals = "/home/stpuser/Results/" + sample_id + "/target.intervals"

	#bait intervals
	cmd1 = "java -jar /software/picard-tools-2.5.0/picard.jar BedToIntervalList I=" + broad_bed + " O=" + bait + " SD=" + "/reference_files/ucsc.hg19.nohap.masked.dict"
	#target intervals
	cmd2 = "java -jar /software/picard-tools-2.5.0/picard.jar BedToIntervalList I=" + small_bed + " O=" + intervals + " SD=" + "/reference_files/ucsc.hg19.nohap.masked.dict"
	print(cmd1)
	print(cmd2)
	run_cmd(cmd1)
	run_cmd(cmd2)
	#generate statistics
	cmd3 = "java -jar /software/picard-tools-2.5.0/picard.jar CollectHsMetrics I=" + in_bam + " O=/home/stpuser/Results/" + sample_id + "/" + sample_id + "hs_metrics.txt R=" + reference + " BAIT_INTERVALS=" + bait + " TARGET_INTERVALS=" + intervals
	print(cmd3)
	run_cmd(cmd3)

	
def call_variants(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_merged_sorted_dupsmarked_RG_indelsrealigned.bam" 
	out_vcf = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_Variants.vcf"
	
	cmd="java -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + reference + " -dt NONE -I " + in_bam + " -stand_call_conf 30 -o " + out_vcf + " -L /home/stpuser/Results/" + sample_id + "\/target.intervals"
	print(cmd)
	run_cmd(cmd)

def decompose_vcf(sample_id):
	in_vcf = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_Variants.vcf"
	out_vcf = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_Variants_decomposed.vcf"
	
	cmd = "vt decompose " + in_vcf + " -o " + out_vcf
	run_cmd(cmd)
	
def annotate_variants(sample_id):
	in_vcf = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_Variants_decomposed.vcf"
	out_vcf = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_Variants_exac.vcf"
	
	cmd= "java -Xmx4g -jar /software/snpEff/SnpSift.jar annotate -v /reference_files/exac/ExAC.r0.3.1.sites.vep.vcf.gz " + in_vcf + " > " + out_vcf
	print(cmd)
	run_cmd(cmd)
	
def main():
	#run_fastqc(sample_list)
	#run_alignment(sample_id, sample_list)
	#qc_read_mappings(sample_id)
	#mark_dups(sample_id)
	#add_readgroup(sample_id)
	#realign_indels(sample_id)
	#get_hs_stats(sample_id) 
	#find_coverage(sample_id)
	#call_variants(sample_id)
	#decompose_vcf(sample_id)
	annotate_variants(sample_id)

if __name__ == '__main__':
    main()
