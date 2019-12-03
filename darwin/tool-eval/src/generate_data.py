#! /usr/bin/python

import os;
import sys;
import glob;
import math;
import subprocess;
import fastqparser;
import numpy as np;

from basicdefines import *;
import fastqparser;



def peek(fp, num_chars):
	data = fp.read(num_chars);
	if len(data) == 0:
		return '';
	fp.seek(num_chars * -1, 1);
	return data;

def get_single_read(fp):
	lines = '';
	
	line = fp.readline();
	header = line;
	lines += line;
	next_char = peek(fp, 1);
	
	num_lines = 1;
	while len(next_char) > 0 and next_char != lines[0] or (next_char == '@' and num_lines < 4):
		line = fp.readline();
		lines += line;
		next_char = peek(fp, 1);
		num_lines += 1;
		
	return [header.rstrip(), lines.rstrip()];

def interleave(reads1_path, reads2_path, out_path):
	fp1 = open(reads1_path, 'r');
	fp2 = open(reads2_path, 'r');
	fp_out = open(out_path, 'w');

	num_read_pairs = 0;
	while True:
		[header1, read1] = get_single_read(fp1);
		[header2, read2] = get_single_read(fp2);
		
		if (len(read1) == 0 and len(read2) > 0) or (len(read1) > 0 and len(read2) == 0) or (header1[0:-3] != header2[0:-3]):
			sys.stderr.write(('ERROR: Reads mismatch! Wrong input files, or reads not in correct order! Read #%d!' % num_read_pairs) + '\n');
			sys.stderr.write(('Header 1: "%s"' % header1) + '\n');
			sys.stderr.write(('Header 2: "%s"' % header2) + '\n');
			break;
		if (len(read1) == 0 and len(read2) == 0):
			break;
		
		fp_out.write(read1 + '\n');
		fp_out.write(read2 + '\n');
		
		num_read_pairs += 1;
	
	fp1.close();
	fp2.close();
	fp_out.close();

def CreateFolders(machine_name, genome_filename):
	if not os.path.exists(machine_name + '/' + genome_filename):
		sys.stderr.write(('Creating %s output folders.' % machine_name) + '\n');
		os.makedirs(machine_name + '/' + genome_filename);

def EstimateCoverageForNumReads(genome_path, genome_filename, mean_read_length, num_reads):
	fp_in = None;
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'

	try:
		fp_in = open(complete_genome_path, 'r');
	except IOError:
		sys.stderr.write(('ERROR: Could not open file "%s" for reading!' % complete_genome_path) + '\n');
		exit(1);
	
	total_genome_length = 0;
	
	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		total_genome_length += len(seq);
		
	fp_in.close();
	
	coverage = int(math.ceil((float(num_reads) / (float(total_genome_length) / float(mean_read_length))))) * 3;
	
	#sys.stderr.write((num_reads) + '\n');
	#sys.stderr.write((total_genome_length) + '\n');
	#sys.stderr.write((mean_read_length) + '\n');
	#print (float(total_genome_length) / float(mean_read_length))
	#print float(num_reads)
	#print coverage
	
	return coverage;

def subsample_alignments_from_fastq(reads_path_prefix, subsampled_set):
	reads_path = reads_path_prefix + '.fq';
	reads_path_fasta = reads_path_prefix + '.fa';
	if (os.path.exists(reads_path) == True):
		complete_reads_path = reads_path_prefix + '-complete_dataset.fq';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		sys.stderr.write(('Renaming file "%s" to "%s"...' % (reads_path, complete_reads_path)) + '\n');
		
		os.rename(reads_path, complete_reads_path);
		
		fp_in = open(complete_reads_path, 'r');
		fp_out = open(reads_path, 'w');

		subsampled_headers = [];
		current_subsample = 0;
		num_read_pairs = 0;
		i = 0;
		while (True):
			[header, read] = get_single_read(fp_in);
			
			if (len(read) == 0):
				break;

			if (i == subsampled_set[current_subsample]):
				fp_out.write(read + '\n');
				subsampled_headers.append(header[1:]);
				current_subsample += 1;
				if (current_subsample >= len(subsampled_set)):
					break;
			
			i += 1;
		
		fp_in.close();
		fp_out.close();

		shell_command = 'rm %s' % (complete_reads_path);
		sys.stderr.write('Removing intermediate file: "%s"\n' % complete_reads_path);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		sys.stderr.write('Converting the FASTQ file to FASTA...');
		fastqparser.convert_to_fasta(reads_path, reads_path_fasta);
		sys.stderr.write('done!\n\n');

		return subsampled_headers;
	else:
		sys.stderr.write('ERROR: Reads file "%s" does not exist!' % (reads_path) + '\n');
		exit(1);

def subsample_alignments_from_sam(reads_path_prefix, subsampled_set):
	sam_path = reads_path_prefix + '.sam';
	if (os.path.exists(sam_path) == True):
		complete_sam_path = reads_path_prefix + '-complete_dataset.sam';
		#complete_reads_path = os.path.dirname(reads_path) + '/' + os.path.splitext(os.path.basename(reads_path))[0] + '-complete_dataset' + os.path.splitext(reads_path)[1];
		
		sys.stderr.write(('Renaming file "%s" to "%s"...' % (sam_path, complete_sam_path)) + '\n');
		
		os.rename(sam_path, complete_sam_path);
		
		fp_in = open(complete_sam_path, 'r');
		fp_out = open(sam_path, 'w');
		
		current_subsample = 0;
		current_num_alignments = 0;
		
		sys.stderr.write(('num_alignments_to_extract = %d' % len(subsampled_set)) + '\n');

		subsampled_qnames = [];

		for line in fp_in:
			if (len(line.strip()) == 0 or line.startswith('@') == True):
				fp_out.write(line);
			else:
				if (current_num_alignments == subsampled_set[current_subsample]):
					fp_out.write(line);
					subsampled_qnames.append(line.split('\t')[0]);
					current_subsample += 1;
					if (current_subsample >= len(subsampled_set)):
						break;

				current_num_alignments += 1;
		
		sys.stderr.write(('current_subsample = %d, subsampled_set[current_subsample] = %d' % (current_subsample, subsampled_set[-1]) ) + '\n');

		fp_in.close();
		fp_out.close();

		shell_command = 'rm %s' % (complete_sam_path);
		sys.stderr.write('Removing intermediate file: "%s"\n' % complete_sam_path);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		return subsampled_qnames;
	else:
		sys.stderr.write(('ERROR: Reads file "%s" does not exist!' % (reads_path)) + '\n');
		exit(1);

def subsample_generated_reads(out_file_prefix, num_reads_to_generate):
	num_generated_reads = 0;
	reads_path = out_file_prefix + '.fq';
	if (os.path.exists(reads_path) == True):
		fp = open(reads_path, 'r');
		while (True):
			[header, read] = fastqparser.get_single_read(fp);
			if (len(read) == 0):
				break;
			num_generated_reads += 1;

		sys.stderr.write('num_generated_reads = %d\n' % num_generated_reads);
		sys.stderr.write('num_reads_to_generate = %d\n' % num_reads_to_generate);
		subsampled_set = sorted(np.random.choice(num_generated_reads, num_reads_to_generate, replace=False));

		subsampled_headers = subsample_alignments_from_fastq(out_file_prefix, subsampled_set);
		subsampled_qnames = subsample_alignments_from_sam(out_file_prefix, subsampled_set);

		# Check the validity of subsampled sets (i.e. if the same reads are included).
		# subsampled_headers = sorted(subsampled_headers);
		# subsampled_qnames = sorted(subsampled_qnames);
		sys.stderr.write('Performing sanity checks on subsampled output...\n');
		if (len(subsampled_headers) != len(subsampled_qnames)):
			sys.stderr.write('ERROR: Subsampling of generated reads and generated SAM file did not produce the same sequences in output! Number of reads in these files differs!\n');
			sys.stderr.write('len(subsampled_headers) = %d\n' % len(subsampled_headers));
			sys.stderr.write('len(subsampled_qnames) = %d\n' % len(subsampled_qnames));
			exit(1);

		current_header = 0;
		current_qname = 0;
		while (current_header < len(subsampled_headers) and current_qname < len(subsampled_qnames)):
			if (subsampled_headers[current_header] != subsampled_qnames[current_qname]):
				sys.stderr.write('ERROR: Subsampling of generated reads and generated SAM file did not produce the same sequences in output! Check the ordering of reads in these files!\n');
				sys.stderr.write('subsampled_headers[current_header] = "%s"\n' % subsampled_headers[current_header]);
				sys.stderr.write('subsampled_qnames[current_qname] = "%s"\n' % subsampled_qnames[current_qname]);
				exit(1);
			current_header += 1;
			current_qname += 1;

		sys.stderr.write('All sanity checks passed!\n');
	else:
		sys.stderr.write('ERROR: Reads not generated! Cannot subsample!\n');
		exit(1);

def Generate454(genome_path, genome_filename, fold_coverage=10, mean_fragsize=1500, std_fragsize=10, machine_name='Roche454', num_reads_to_generate=-1):
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	simulator_bin = TOOLS_ROOT_ABS + '/art_bin_VanillaIceCream/art_454';
#	out_file_prefix = machine_name + '/' + genome_filename + '/' + genome_filename;

	is_paired_end = False;
	
	if mean_fragsize==0 or std_fragsize==0:
		machine_name += '-single_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' -s -r 1403002416 ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage);
		shell_command = simulator_bin + r' -s ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage);
		is_paired_end = False;
	else:
		machine_name += '-paired_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' -s -r 1403002416 ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage) + r' ' + str(mean_fragsize) + r' ' + str(std_fragsize);
		shell_command = simulator_bin + r' -s ' + complete_genome_path + r' ' + out_file_prefix + r' ' + str(fold_coverage) + r' ' + str(mean_fragsize) + r' ' + str(std_fragsize);
		is_paired_end = True;
	
	sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');	
	subprocess.call(shell_command, shell=True);

	# if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# 	ExtractNReadsFromFile(out_file_prefix, NUM_READS_TO_GENERATE);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, NUM_READS_TO_GENERATE);
	# if (num_reads_to_generate > 0):
	# 	ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	if (num_reads_to_generate > 0):
		subsample_generated_reads(out_file_prefix, num_reads_to_generate);
	
	if is_paired_end == True:
		sys.stderr.write(('Interleaving paired end reads to file "%s"' % (out_file_prefix + '.fq')) + '\n');
		interleave(out_file_prefix + '1.fq', out_file_prefix + '2.fq', out_file_prefix + '.fq');

def GenerateIllumina(genome_path, genome_filename, read_length=100, fold_coverage=10, mean_fragsize=500, std_fragsize=10, machine_name='Illumina', num_reads_to_generate=-1):
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	simulator_bin = TOOLS_ROOT_ABS + '/art_bin_VanillaIceCream/art_illumina';
	
	is_paired_end = False;
	
	if mean_fragsize==0 or std_fragsize==0:			# Single end reads
		machine_name += '-single_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' --rndSeed 1403000281 -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -sam';	# Single-end
		shell_command = simulator_bin + r' -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -sam';	# Single-end
		is_paired_end = False;
	else:			# Paired end reads
		machine_name += '-paired_end';
		CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
		CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
		
		out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
		#shell_command = simulator_bin + r' --rndSeed 1403000281 -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -m ' + str(mean_fragsize) + r' -s ' + str(std_fragsize) + r' -sam';	# Paired-end
		shell_command = simulator_bin + r' -i ' + complete_genome_path + r' -o ' + out_file_prefix + r' -l ' + str(read_length) + r' -f ' + str(fold_coverage) + r' -m ' + str(mean_fragsize) + r' -s ' + str(std_fragsize) + r' -sam';	# Paired-end
		is_paired_end = True;

	sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
	subprocess.call(shell_command, shell=True);
	
	# if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# 	ExtractNReadsFromFile(out_file_prefix, NUM_READS_TO_GENERATE);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, NUM_READS_TO_GENERATE);
	# if (num_reads_to_generate > 0):
	# 	ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	if (num_reads_to_generate > 0):
		subsample_generated_reads(out_file_prefix, num_reads_to_generate);
	
	if is_paired_end == True:
		sys.stderr.write(('Interleaving paired end reads to file "%s"' % (out_file_prefix + '.fq')) + '\n');
		interleave(out_file_prefix + '1.fq', out_file_prefix + '2.fq', out_file_prefix + '.fq');
	


def GetPBSimRefName(ref_file):
	try:
		fp = open(ref_file, 'r');
	except IOError:
		sys.stderr.write(('ERROR: Could not open file "%s" for reading!' % ref_file) + '\n');
		exit(1);

	header = fp.readline();
	header = header.rstrip();
	fp.close();
	
	if not header.startswith('>'):
		sys.stderr.write(("ERROR: PBsim's ref file does not start with a FASTA header!") + '\n');
		exit(1);

	ref_name = header[1:];
	trimmed_ref_name = ref_name.split()[0];
	
	return [ref_name, trimmed_ref_name];

# --difference-ratio   ratio of differences. substitution:insertion:deletion.
def GeneratePacBio(genome_path, genome_filename, fold_coverage=30, length_mean=10000, length_sd=0, length_min=10000, length_max=10000, accuracy_mean=0.78, accuracy_sd=0.02, accuracy_min=0.75, difference_ratio='10:60:30', machine_name='PacBio', num_reads_to_generate=-1):
#	machine_name = 'PacBio';
	CreateFolders(READS_SIMULATED_ROOT_ABS + '/' + machine_name, genome_filename);
	CreateFolders(EVALUATION_PATH_ROOT_ABS + '/' + READS_SIMULATED_ROOT + '/' + machine_name, genome_filename);
	
	complete_genome_path = genome_path + '/' + genome_filename + '.fa'
	out_file_prefix = READS_SIMULATED_ROOT_ABS + '/' + machine_name + '/' + genome_filename + '/' + 'reads';
	
	simulator_path = TOOLS_ROOT_ABS + '/pbsim-1.0.3-Linux-amd64';
	simulator_bin = TOOLS_ROOT_ABS + '/pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim';
	
	final_sam_file = out_file_prefix + '.sam';
	fp = open(final_sam_file, 'w');
	fp.close();
	
	final_fastq_file = out_file_prefix + '.fq';
	fp = open(final_fastq_file, 'w');
	fp.close();
# SIMULATOR_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/tools/pbsim-1.0.3-Linux-amd64
# REFERENCE_PATH=/home/ivan/work/eclipse-workspace/golden-bundle/reference-genomes/escherichia_coli.fa
# $SIMULATOR_PATH/Linux-amd64/bin/pbsim --data-type CLR --depth 20 --model_qc $SIMULATOR_PATH/data/model_qc_clr --seed 1234567890 --prefix $OUT_FILE_PREFIX $REFERENCE_PATH
# /last-460/scripts/maf-convert.py sam $MAF_FILE $SAM_FILE

	#random_seed = '1234567890';
	random_seed = '32874638';

	# Data type:
	#	Continuous Long Read (CLR) : long and high error rate.
	#	Circular consensus Read (CCS) : short and low error rate.
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	shell_command = simulator_bin + r' --data-type CLR --depth ' + str(fold_coverage) + \
					' --model_qc ' + simulator_path + \
					'/data/model_qc_clr ' + \
					' --length-mean ' + str(length_mean) + \
					' --length-sd ' + str(length_sd) + \
					' --length-min ' + str(length_min) + \
					' --length-max ' + str(length_max) + \
					' --accuracy-mean ' + str(accuracy_mean) + \
					' --accuracy-sd ' + str(accuracy_sd) + \
					' --accuracy-min ' + str(accuracy_min) + \
					' --difference-ratio ' + difference_ratio + \
					' --prefix ' + out_file_prefix + \
					r' ' + complete_genome_path;
				
					#' --seed ' + random_seed + \
	
	sys.stderr.write(('Simulating PacBio reads using PBsim') + '\n');
	sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
	
	#exit(1);
	subprocess.call(shell_command, shell=True);

	sys.stderr.write((' ') + '\n');
	
	sys.stderr.write(('Converting generated *.maf files to corresponding SAM files.') + '\n');

	maf_files = glob.glob(out_file_prefix + '*.maf');
	maf_files = sorted(maf_files);
	sys.stderr.write('\n'.join(maf_files) + '\n');
	
	for maf_file in maf_files:		
		# Convert the maf file to SAM format.
		sam_file = maf_file[0:-3] + 'sam';
		shell_command = MAF_CONVERT_ROOT_ABS + '/maf-convert sam ' + maf_file + ' > ' + sam_file;
		sys.stderr.write(('Converting MAF to SAM ("%s" -> "%s")' % (maf_file, sam_file)) + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		
		# Use Bash's sed to replace the 'ref' keyword in the generated maf files with the actual name of the reference sequence, and concatenate all the separate SAM files to one file.
		[reference_name, trimmed_reference_name] = GetPBSimRefName(maf_file[0:-3] + 'ref');
		# Here we escape the special characters so that SED command runs properly if any of these characters should to appear in a FASTA header.
		escape_chars = r'\/()[].*^$';
		reference_name = ''.join([('\\' + char) if char in escape_chars else char for char in reference_name]);
		shell_command = r'cat ' + sam_file + r" | sed 's/^\(.*\)ref/\1" + reference_name + r"/' >> " + final_sam_file
		sys.stderr.write(('Replacing PBsim\'s "ref" keyword with actual FASTA header') + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		
		fastq_file = maf_file[0:-3] + 'fastq';
		shell_command = r'cat ' + fastq_file + ' >> ' + final_fastq_file;
		sys.stderr.write(('Concatenating FASTQ file to the total reads file ("%s" -> "%s")' % (fastq_file, final_fastq_file)) + '\n');
		sys.stderr.write(('Executing command: "%s"' % shell_command) + '\n');
		subprocess.call(shell_command, shell=True);
		sys.stderr.write((' ') + '\n');
		
		sys.stderr.write((' ') + '\n');

		shell_command = 'rm %s' % (maf_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % maf_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);
		shell_command = 'rm %s' % (sam_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % sam_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);
		shell_command = 'rm %s' % (fastq_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % fastq_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

		ref_file = maf_file[0:-3] + 'ref';
		shell_command = 'rm %s' % (ref_file);
		sys.stderr.write('Removing intermediate file: "%s"\n' % ref_file);
		sys.stderr.write(('Executing command: "%s"\n\n' % shell_command));
		subprocess.call(shell_command, shell=True);

#	if (GENERATE_FIXED_AMOUNT_OF_READS == True):
	# if (num_reads_to_generate > 0):
	# 	ExtractNReadsFromFile(out_file_prefix, num_reads_to_generate);
	# 	ExtractNAlignmentsFromSAM(out_file_prefix, num_reads_to_generate);
	if (num_reads_to_generate > 0):
		subsample_generated_reads(out_file_prefix, num_reads_to_generate);

	#sam_files = glob.glob(out_file_prefix + '*.sam');
	#sam_files = sorted(sam_files);
	#sys.stderr.write((sam_files) + '\n');

def GenerateOxfordNanoporeFromObservedStatistics(genome_filename, num_reads_to_generate=-1):
	if (num_reads_to_generate <= 0):
		coverage = 30;
		machine_suffix = '-cov30';
	else:
		# num_reads_to_generate = 10000;
		mean_read_length = 1000;
		coverage = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, genome_filename, mean_read_length, num_reads_to_generate) + 1;
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);

	## The maximum value for length_max parameter (limited by PBsim) is 100000.

	sys.stderr.write(('num_reads_to_generate = %d' % num_reads_to_generate) + '\n');

	# These are simulations for difference_ratio obtained from LAST's alignments of E. Coli R7.3 reads (from Nick Loman).
	# 1d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.41	0.05	0.40	0.10	0.60
	# Insertion rate stats: 	0.05	0.02	0.05	0.00	0.23
	# Deletion rate stats:  	0.16	0.05	0.15	0.00	0.49
	# Mismatch rate stats:  	0.20	0.03	0.20	0.03	0.32
	# Match rate stats:     	0.75	0.04	0.75	0.59	0.97
	# Read length stats:    	3629.76	3294.04	2438.00	57.00	31299.00
	# Difference ratio: 51:11:38 (mismatch:insertion:deletion)
	#
	# This simulates data with realistic observed error rate:
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=100000, accuracy_mean=(1.0 - 0.41), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38', machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# ISCA GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=10000, length_sd=0, length_min=10000, length_max=10000, accuracy_mean=(1.0 - 0.41), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38', machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# The following call would simulate a generic 40% error rate (not realistic): # GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.40), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38',
	# 															machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# 2d reads:
	# [CIGAR statistics - individual indels]
	#                       	mean	std	median	min	max
	# Error rate stats:     	0.31	0.09	0.30	0.02	0.59
	# Insertion rate stats: 	0.05	0.03	0.05	0.00	0.28
	# Deletion rate stats:  	0.09	0.06	0.08	0.00	0.48
	# Mismatch rate stats:  	0.16	0.07	0.16	0.00	0.36
	# Match rate stats:     	0.78	0.07	0.79	0.59	0.99
	# Read length stats:    	2006.14	3015.25	614.00	42.00	28601.00
	# Difference ratio: 55:17:28 (mismatch:insertion:deletion)
	#
	# This simulates data with realistic observed error rate:
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=5600, length_sd=3500, length_min=100, length_max=100000, accuracy_mean=(1.0 - 0.31), accuracy_sd=0.09, accuracy_min=(1.0 - 0.60), difference_ratio='55:17:28', machine_name='OxfordNanopore-pbsim-observed_last-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=10000, length_sd=0, length_min=10000, length_max=10000, accuracy_mean=(1.0 - 0.31), accuracy_sd=0.09, accuracy_min=(1.0 - 0.60), difference_ratio='55:17:28', machine_name='OxfordNanopore-pbsim-observed_last-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# The following call would simulate a generic 20% error rate (not realistic):
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.20), accuracy_sd=0.05, accuracy_min=(1.0 - 0.50), difference_ratio='55:17:28',
	# 															machine_name='OxfordNanopore-pbsim-observed_last-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);



	# These are simulations for difference_ratio obtained from GraphMap's alignments of E. Coli R7.3 reads (from Nick Loman).
	# 1d reads:
	# [CIGAR statistics - individual indels] Newest (from 20150329).
	#                         mean    std     median  min     max
	# Error rate:             0.44    0.05    0.46    0.25    0.97
	# Insertion rate:         0.15    0.07    0.14    0.00    0.38
	# Deletion rate:          0.11    0.06    0.11    0.00    0.28
	# Mismatch rate:          0.18    0.03    0.18    0.05    0.71
	# Match rate:             0.67    0.10    0.68    0.17    0.88
	# Read length:            4713.15 3937.17 4026.00 80.00   94116.00
	# Difference ratio: 42:34:24 (mismatch:insertion:deletion)
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=10000, length_sd=0, length_min=10000, length_max=10000, accuracy_mean=(1.0 - 0.44), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='42:34:24', 	machine_name='OxfordNanopore-pbsim-observed_graphmap-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate); 
        # The following call would simulate a generic 40% error rate (not realistic):
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.40), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='44:29:27',
	# 															machine_name='OxfordNanopore-pbsim-observed_graphmap-1d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);

	# 2d reads:
	# [CIGAR statistics - individual indels] Newest (from 20150329).
	#                         mean    std     median  min     max
	# Error rate:             0.32    0.11    0.28    0.13    0.90
	# Insertion rate:         0.11    0.07    0.07    0.03    0.34
	# Deletion rate:          0.10    0.05    0.08    0.00    0.33
	# Mismatch rate:          0.11    0.06    0.09    0.03    0.68
	# Match rate:             0.78    0.12    0.84    0.24    0.93
	# Read length:            5651.16 3466.81 5316.00 102.00  33550.00
	# Difference ratio: 37:33:30 (mismatch:insertion:deletion)
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=5600, length_sd=3500, length_min=100, length_max=100000,
	# 															accuracy_mean=(1.0 - 0.32), accuracy_sd=0.11, accuracy_min=(1.0 - 0.60), difference_ratio='37:33:30',
	# 															machine_name='OxfordNanopore-pbsim-observed_graphmap-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	# The following call would simulate a generic 20% error rate (not realistic):
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=94000,
	# 															accuracy_mean=(1.0 - 0.20), accuracy_sd=0.05, accuracy_min=(1.0 - 0.50), difference_ratio='37:23:40',
	# 															machine_name='OxfordNanopore-pbsim-observed_graphmap-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);

def GenerateOxfordNanopore2dFromReportedStatistics(genome_filename, num_reads_to_generate=-1):
	# Paper reference:
	# Miten Jain, Ian T. Fiddes, Karen H Miga, Hugh E. Olsen, Benedict Paten & Mark Akeson:
	# Improved data analysis for the MinION nanopore sequencer. Nature Methods (2015) doi:10.1038/nmeth.3290
	# http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3290.html

	# "We evaluated and optimized the performance of the MinION nanopore sequencer using M13 genomic DNA and used
	# expectation maximization to obtain robust maximum-likelihood estimates for insertion, deletion and substitution
	# error rates (4.9%, 7.8% and 5.1%, respectively). Over 99% of high-quality 2D MinION reads mapped to the
	# reference at a mean identity of 85%. 
	# ...
	# "This showed that insertions were less frequent than deletions by about twofold in 2D reads and about threefold in template and complement reads. The combined
	# insertion-deletion (indel) rate was between 0.13 (2D reads) and 0.2 (template and complement reads) events per aligned base. For all
	# read types, indels were predominantly single bases (Supplementary Fig. 6). Substitutions varied from 0.21 (for template reads) to 0.05
	# (for 2D reads) events per aligned base (Fig. 3c and Supplementary Figs. 7 and 8). Substitution errors were not uniform; in particular,
	# A-to-T and T-to-A errors were estimated to be very low, at 0.04% and 0.1%, respectively (Supplementary Note 1)."

	error_rate = (0.049 + 0.078 + 0.051);	# = 0.178
	insertion_ratio = int((0.049 / error_rate) * 100);
	deletion_ratio = int((0.078 / error_rate) * 100);
	mismatch_ratio = 100 - insertion_ratio - deletion_ratio;
	difference_ratio = '%d:%d:%d' % (mismatch_ratio, insertion_ratio, deletion_ratio);	

	if (num_reads_to_generate <= 0):
		coverage = 20;
		machine_suffix = '-cov20';
	else:
		mean_read_length = 6200;
		coverage = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, genome_filename, mean_read_length, num_reads_to_generate) + 1;
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);

	sys.stderr.write(('num_reads_to_generate = %d' % num_reads_to_generate) + '\n');

	## The maximum value for length_max parameter (limited by PBsim) is 100000.
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=5600, length_sd=3500, length_min=100, length_max=100000,
																accuracy_mean=(1.0 - (0.049 + 0.078 + 0.051)), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio=difference_ratio,
																machine_name='OxfordNanopore-pbsim-observed_marginalign-2d' + machine_suffix, num_reads_to_generate=num_reads_to_generate);



def GenerateNGSData(num_reads_to_generate=-1):
	##### NGS DATA ###
	if (num_reads_to_generate <= 0):
		coverage_nmeni = 20;
		coverage_ecoli = 20;
		coverage_scerevisiae = 20;
		coverage_celegans = 20;
		coverage_hg19v38chr3 = 20;
		machine_suffix = '-cov20';
	else:
		mean_read_length = 150;
		coverage_nmeni = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for neisseria_meningitidis: %d' % coverage_nmeni) + '\n');
		coverage_ecoli = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for escherichia_coli: %d' % coverage_ecoli) + '\n');
		coverage_scerevisiae = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for saccharomyces_cerevisiae: %d' % coverage_scerevisiae) + '\n');
		coverage_celegans = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for caenorhabditis_elegans: %d' % coverage_celegans) + '\n');
		coverage_hg19v38chr3 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for hg19_v38-chr3: %d' % coverage_hg19v38chr3) + '\n');
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);
		coverage_hg19v38 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for hg19: %d' % coverage_hg19v38) + '\n');
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);

	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', read_length=150, fold_coverage=coverage_nmeni, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', read_length=150, fold_coverage=coverage_ecoli, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', read_length=150, fold_coverage=coverage_scerevisiae, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', read_length=150, fold_coverage=coverage_celegans, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', read_length=150, fold_coverage=coverage_hg19v38chr3, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GenerateIllumina(REFERENCE_GENOMES_ROOT_ABS, 'hg19', read_length=150, fold_coverage=coverage_hg19v38, mean_fragsize=0, std_fragsize=0, machine_name='Illumina' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', fold_coverage=coverage_nmeni, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=coverage_ecoli, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', fold_coverage=coverage_scerevisiae, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', fold_coverage=coverage_celegans, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	#Generate454(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', fold_coverage=coverage_hg19v38chr3, mean_fragsize=0, std_fragsize=0, machine_name='Roche454' + machine_sufix);
	##################

def GeneratePacBioData(num_reads_to_generate=-1):
	##### PACBIO NORMAL DATA #####
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=80, length_mean=3000, length_sd=2300.0, accuracy_mean=0.78, accuracy_sd=0.02, accuracy_min=0.75, machine_name='PacBio');
	if (num_reads_to_generate <= 0):
		coverage_nmeni = 20;
		coverage_ecoli = 30;
		coverage_scerevisiae = 20;
		coverage_celegans = 30;
		coverage_hg19v38chr3 = 20;
		machine_suffix = '-cov30';
	else:
		mean_read_length = 1000;
		coverage_nmeni = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', mean_read_length, num_reads_to_generate) + 1;
		#sys.stderr.write(('Coverage for neisseria_meningitidis: %d' % coverage_nmeni) + '\n');
		coverage_ecoli = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for escherichia_coli: %d' % coverage_ecoli) + '\n');
		coverage_scerevisiae = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', mean_read_length, num_reads_to_generate) + 1;
		#sys.stderr.write(('Coverage for saccharomyces_cerevisiae: %d' % coverage_scerevisiae) + '\n');
		coverage_celegans = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', mean_read_length, num_reads_to_generate) + 1;
		sys.stderr.write(('Coverage for caenorhabditis_elegans: %d' % coverage_celegans) + '\n');
		coverage_hg19v38chr3 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', mean_read_length, num_reads_to_generate) + 1;
		#sys.stderr.write(('Coverage for hg19_v38-chr3: %d' % coverage_hg19v38chr3) + '\n');
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);
		coverage_hg19v38 = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38', mean_read_length, num_reads_to_generate) + 1;
		#sys.stderr.write(('Coverage for hg19: %d' % coverage_hg19v38) + '\n');
		machine_suffix = '-%dk' % (num_reads_to_generate / 1000);
		
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'neisseria_meningitidis', fold_coverage=coverage_nmeni, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=coverage_ecoli, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'saccharomyces_cerevisiae', fold_coverage=coverage_scerevisiae, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'caenorhabditis_elegans', fold_coverage=coverage_celegans, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'hg19_v38-chr3', fold_coverage=coverage_hg19v38chr3, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	#GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'hg19', fold_coverage=coverage_hg19v38, machine_name='PacBio' + machine_suffix, num_reads_to_generate=num_reads_to_generate);
	##############################

def GenerateOxfordNanoporeDataObserved(num_reads_to_generate=-1):
	##### OXFORD NANOPORE DATA #####
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.

	#GenerateOxfordNanoporeFromObservedStatistics('neisseria_meningitidis', num_reads_to_generate=num_reads_to_generate);
	#GenerateOxfordNanoporeFromObservedStatistics('escherichia_coli', num_reads_to_generate=num_reads_to_generate);
	#GenerateOxfordNanoporeFromObservedStatistics('saccharomyces_cerevisiae', num_reads_to_generate=num_reads_to_generate);
	#GenerateOxfordNanoporeFromObservedStatistics('caenorhabditis_elegans', num_reads_to_generate=num_reads_to_generate);
	#GenerateOxfordNanoporeFromObservedStatistics('hg19_v38-chr3', num_reads_to_generate=num_reads_to_generate);
	#GenerateOxfordNanoporeFromObservedStatistics('hg19', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanoporeFromObservedStatistics('human_g1k_v37', num_reads_to_generate=num_reads_to_generate);

def GenerateOxfordNanoporeDataMarginAlign(num_reads_to_generate=-1):
	##### OXFORD NANOPORE DATA #####
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.	
	GenerateOxfordNanopore2dFromReportedStatistics('neisseria_meningitidis', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanopore2dFromReportedStatistics('escherichia_coli', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanopore2dFromReportedStatistics('saccharomyces_cerevisiae', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanopore2dFromReportedStatistics('caenorhabditis_elegans', num_reads_to_generate=num_reads_to_generate);
	GenerateOxfordNanopore2dFromReportedStatistics('hg19_v38-chr3', num_reads_to_generate=num_reads_to_generate);

def GenerateGridTest(num_reads_to_generate):
	##### OXFORD NANOPORE DATA #####
	# --difference-ratio   ratio of differences. substitution:insertion:deletion.
	# GenerateOxfordNanoporeFromObservedStatistics('caenorhabditis_elegans', num_reads_to_generate=num_reads_to_generate);
	
	# genome_filename = 'caenorhabditis_elegans';
	genome_filename = 'saccharomyces_cerevisiae';
	error_rates = [0.05, 0.10, 0.15, 0.20, 0.25];
	read_lengths = [1000, 2000, 3000, 4000, 5000];
	# error_rates = [0.05, 0.10];
	# read_lengths = [1000, 2000];
#	error_rates = [0.05];
#	read_lengths = [1000];

	sys.stderr.write(('num_reads_to_generate = %d' % num_reads_to_generate) + '\n');

	for read_length in read_lengths:
		for error_rate in error_rates:
			mean_read_length = read_length;
			coverage = EstimateCoverageForNumReads(REFERENCE_GENOMES_ROOT_ABS, genome_filename, mean_read_length, num_reads_to_generate) + 1;
			machine_suffix = '-l%dk-e%.2f' % ((read_length / 1000), error_rate);

			GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, genome_filename, fold_coverage=coverage, length_mean=read_length, length_sd=1, length_min=50, length_max=100000,
																		accuracy_mean=(1.0 - error_rate), accuracy_sd=0.01, accuracy_min=(1.0 - error_rate - 0.01), difference_ratio='51:11:38',
																		machine_name='GridSearch' + machine_suffix, num_reads_to_generate=num_reads_to_generate);


def GenerateAll():
	# num_reads_to_generate = 10000;
	#num_reads_to_generate = 1000;

	#GenerateNGSData(num_reads_to_generate);
	#GeneratePacBioData(num_reads_to_generate);
	#GenerateOxfordNanoporeDataObserved(num_reads_to_generate);
	GenerateOxfordNanoporeDataObserved();
	#GeneratePacBioData();

	#GenerateGridTest(10000);


	# num_reads_to_generate = 1;
	# GenerateOxfordNanoporeDataMarginAlign(num_reads_to_generate);
	# GenerateNGSData(num_reads_to_generate);

	# coverage = 60;
	# machine_suffix = '-cov60';
	# GeneratePacBio(REFERENCE_GENOMES_ROOT_ABS, 'escherichia_coli', fold_coverage=coverage, length_mean=4400, length_sd=3900, length_min=50, length_max=100000,
	# 															accuracy_mean=(1.0 - 0.41), accuracy_sd=0.05, accuracy_min=(1.0 - 0.60), difference_ratio='51:11:38',
	# 															machine_name='OxfordNanopore-pbsim-observed_last-1d' + machine_suffix, num_reads_to_generate=-1);



def Main():
	GenerateAll();

if __name__ == '__main__':
	Main();
