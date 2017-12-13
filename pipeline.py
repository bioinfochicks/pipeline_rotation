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
	
def run_alignment(r1, r2, sample_id):
	cmd = "bwa mem /reference_files/ucsc.hg19.nohap.masked.fasta {0} {1} -v 2 > /home/stpuser/Results/{2}/{2}_aln.sam".format(r1, r2, sample_id)
	run_cmd(cmd)
	return

def convert_sam_bam(sample_id):
	in_sam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_aln.sam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_sorted.bam"
	temp_bam = "/home/stpuser/Results/" + sample_id + "/temp.bam"
	
	cmd = "samtools view -bS " + in_sam + "  | samtools sort - | samtools index - -o " + out_bam
	run_cmd(cmd)
	os.remove(in_sam)
	return 
	
def qc_read_mappings(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_sorted.bam"
	out_text = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_samflagstats.txt"
	cmd1 = "samtools flagstat " + in_bam + " > " + out_text
	run_cmd(cmd1) 
	

	
def mark_dups(sample_id):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_sorted.bam"
	out_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_sorted_dupsmarked.bam"
	cmd = "sambamba markdup -p " + in_bam + " " + out_bam
	run_cmd(cmd)
	
def find_coverage(sample_id, bed):
	in_bam = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_sorted_dupsmarked.bam"
	out1 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_coverage_depth_base_small_panel.txt"
	out2 = "/home/stpuser/Results/" + sample_id + "/" + sample_id + "_coverage_depth_region_small_panel.txt"
	print("Input BAM file: " + in_bam)
	print("...running sambamba depth base...")
	cmd1 = "sambamba depth base --min-coverage=0 -q29 -m  -L " + bed + " " + in_bam + " > " + out1 
	print("Finished sambamba depth base...")
	run_cmd(cmd1) 
	print("...running sambamba region base...")
	cmd2 = "sambamba depth region --cov-threshold=0 --cov-threshold=10 --cov-threshold=20 --cov-threshold=30 --cov-threshold=40 --cov-threshold=50 -q29 -m -L " + bed + " > " + out2
	run_cmd(cmd2)
	print("Finished sambamba depth region...")

	
#def realign_indels(sample_id):

def main():
	#run_fastqc(sample_list)
	#concat_files(sample_id, sample_list, out_dir)
	#r1 = "/home/stpuser/Results/" + sample_id + "/R1.fastq.gz"
	#r2 = "/home/stpuser/Results/" + sample_id + "/R2.fastq.gz"
	#run_alignment(r1, r2, sample_id)
	#convert_sam_bam(sample_id)
	#qc_read_mappings(sample_id)
	#mark_dups(sample_id)
	find_coverage(sample_id, small_bed)

if __name__ == '__main__':
    main()

