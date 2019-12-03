#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

import basicdefines;

ALIGNER_URL = 'https://github.com/PacificBiosciences/blasr.git';
ALIGNER_PATH = SCRIPT_PATH + '/../aligners/blasr';
BIN = 'alignment/bin/blasr';
MAPPER_NAME = 'BLASR';



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#	reads_file			Path to a FASTA/FASTQ file containing reads.
#	reference_file		Path to a reference genome FASTA file.
#	machine_name		A symbolic name to specify a set of parameters for a specific sequencing platform.
#	output_path			Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#	output_suffix		A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
	parameters = '';
	num_threads = multiprocessing.cpu_count() / 2;

	if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
		parameters = '-nproc %s -sam -bestn 1 -minMatch 7' % str(num_threads);

	elif ((machine_name.lower() == 'pacbio')):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);

	elif ((machine_name.lower() == 'nanopore')):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);

	elif ((machine_name.lower() == 'longindel')):
		parameters = '-nproc %s -sam -bestn 1 -clipping none' % str(num_threads);
   # -clipping [none|hard|subread|soft] (none)
   #             Use no/hard/subread/soft clipping for SAM output.

	elif ((machine_name.lower() == 'debug')):
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);

	elif ((machine_name.lower() == 'pacbiom4')):
		parameters = '-nproc %s -bestn 1 -m 4' % str(num_threads);

	else:			# default
		parameters = '-nproc %s -sam -bestn 1' % str(num_threads);



	if (output_suffix != ''):
		output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;
	
	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, output_filename);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);
	
	# Run the indexing process, and measure execution time and memory.
	if ((not os.path.exists(reference_file + '.blasrsa'))):
		sys.stderr.write('[%s wrapper] Generating index...\n' % (MAPPER_NAME));
		command = '%s %s/alignment/bin/sawriter %s.blasrsa %s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_PATH, reference_file, reference_file);
		sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
		subprocess.call(command, shell=True);
		sys.stderr.write('\n\n');
	else:
		sys.stderr.write('[%s wrapper] Reference index already exists. Continuing.\n' % (MAPPER_NAME));
		sys.stderr.flush();

	# Run the alignment process, and measure execution time and memory.
	sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
	command = '%s %s/%s %s %s %s -sa %s.blasrsa -out %s' % (basicdefines.measure_command(memtime_file), ALIGNER_PATH, BIN, reads_file, reference_file, parameters, reference_file, sam_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');
	
	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));
	
	if (not os.path.exists(basicdefines.ALIGNERS_PATH_ROOT_ABS)):
		os.makedirs(basicdefines.ALIGNERS_PATH_ROOT_ABS);

	sys.stderr.write('[%s wrapper] Cloning git repository.\n' % (MAPPER_NAME));
	command = 'cd %s; git clone %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, ALIGNER_URL);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	# sys.stderr.write('[%s wrapper] Checking out commit "0087df7bce410b2f9cb7d03efff649106c800b28" for reproducibility purposes.\n' % (MAPPER_NAME));
	# command = 'cd %s; git checkout 0087df7bce410b2f9cb7d03efff649106c800b28' % (ALIGNER_PATH);
	sys.stderr.write('[%s wrapper] Checking out commit "f7bf1e56871d747829a6a34b13b50debdebf1d0b" for reproducibility purposes.\n' % (MAPPER_NAME));
	command = 'cd %s; git checkout f7bf1e56871d747829a6a34b13b50debdebf1d0b' % (ALIGNER_PATH);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	yes_no = raw_input("[%s wrapper] requires some libraries to be installed. Continue? [y/n] " % (MAPPER_NAME));
	if (yes_no != 'y'):
		return;

	sys.stderr.write('[%s wrapper] Please note that the installation of these libraries assumes that the OS is Ubuntu/Debian based.\n' % (MAPPER_NAME));
	sys.stderr.write('[%s wrapper] Sudo will be required.\n' % (MAPPER_NAME));

	command = 'sudo apt-get install libhdf5-dev';
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	sys.stderr.write('[%s wrapper] Running make.\n' % (MAPPER_NAME));
	command = 'cd %s; make' % (ALIGNER_PATH);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');
	
	sys.stderr.write('[%s wrapper] All instalation steps finished.\n' % (MAPPER_NAME));
	sys.stderr.write('\n');

# The newer version is a bit more complicated to build.
# mkdir deps; cd deps; wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/bin/linux-x86_64-gnu481/hdf5-1.8.12-linux-x86_64-gnu481-shared.tar.gz;
# tar -xvf hdf5-1.8.12-linux-x86_64-gnu481-shared.tar.gz; cd ..;
# ./configure.py --no-pbbam HDF5_INCLUDE=$ROOT_PATH/aligneval/aligners/blasr/deps/hdf5-1.8.12-linux-x86_64-gnu481-shared/include HDF5_LIB=$ROOT_PATH/aligneval/aligners/blasr/deps/hdf5-1.8.12-linux-x86_64-gnu481-shared/lib;
# make build-submodule; make



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\t- mode          - either "run" or "install". If "install" other parameters can be ommitted.\n');
	sys.stderr.write('\t- machine_name  - "illumina", "roche", "pacbio", "nanopore" or "default".\n');
	sys.stderr.write('\t- output_suffix - suffix for the output filename.\n');

	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) < 2 or len(sys.argv) > 7):
		verbose_usage_and_exit();

	if (sys.argv[1] == 'install'):
		download_and_install();
		exit(0);

	elif (sys.argv[1] == 'run'):
		if (len(sys.argv) < 6):
			verbose_usage_and_exit();

		reads_file = sys.argv[2];
		reference_file = sys.argv[3];
		machine_name = sys.argv[4];
		output_path = sys.argv[5];
		output_suffix = '';

		if (len(sys.argv) == 7):
			output_suffix = sys.argv[6];
		run(reads_file, reference_file, machine_name, output_path, output_suffix);

	else:
		verbose_usage_and_exit();
