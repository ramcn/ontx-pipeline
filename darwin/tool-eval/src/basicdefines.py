 #! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;

import fnmatch
import subprocess;



ALIGNERS_PATH = 'aligners';
ALIGNERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + ALIGNERS_PATH;

EVALUATION_PATH = 'evaluation';
EVALUATION_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + EVALUATION_PATH;

RESULTS_PATH = 'results';
RESULTS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + RESULTS_PATH;

WRAPPERS_PATH = 'wrappers';
WRAPPERS_PATH_ROOT_ABS = SCRIPT_PATH + '/../' + WRAPPERS_PATH;

REFERENCE_GENOMES_ROOT = 'reference-genomes';
REFERENCE_GENOMES_ROOT_ABS = SCRIPT_PATH + '/../' + REFERENCE_GENOMES_ROOT;



READS_REAL_ROOT = 'reads-real';
READS_REAL_ROOT_ABS = SCRIPT_PATH + '/../' + READS_REAL_ROOT;

READS_SIMULATED_ROOT = 'reads-simulated';
READS_SIMULATED_ROOT_ABS = SCRIPT_PATH + '/../' + READS_SIMULATED_ROOT;



# INDEX_ROOT = 'index';
# INDEX_ROOT_ABS = SCRIPT_PATH + '/../' + INDEX_ROOT;

# ALIGNMENTS_ROOT = 'alignments_for_testing';
# ALIGNMENTS_ROOT_ABS = SCRIPT_PATH + '/../' + ALIGNMENTS_ROOT;
# ALIGNMENTS_ROOT = 'alignments_for_testing';
# ALIGNMENTS_ROOT_ABS = SCRIPT_PATH + '/../' + ALIGNMENTS_ROOT;

# RESULTS_ROOT = 'results';
# RESULTS_ROOT_ABS = SCRIPT_PATH + '/../' + RESULTS_ROOT;

TOOLS_ROOT = 'tools';
TOOLS_ROOT_ABS = SCRIPT_PATH + '/../' + TOOLS_ROOT;

MAF_CONVERT_ROOT = TOOLS_ROOT + '/last-534/scripts';
MAF_CONVERT_ROOT_ABS = SCRIPT_PATH + '/../' + MAF_CONVERT_ROOT;






# Finds all files with a given filename pattern starting from the given path, recursivelly.
def find_files(start_path, file_pattern, max_depth=-1):
	matches = []
	for root, dirnames, filenames in os.walk(start_path):
		for filename in fnmatch.filter(filenames, file_pattern):
			depth = len(root.split('/'));
			if (max_depth < 0 or (max_depth >= 0 and depth <= max_depth)):
				matches.append(os.path.join(root, filename))

	matches.sort();
	return matches;

# Simply lists all subfolders, recursivelly. Parameter 'depth' can be used to specify the minimum depth of the subfolders.
def find_folders(start_path, depth=0):
	matches = []
	for x in os.walk(start_path):
		if (len(x[0].split('/')) >= depth):
			matches.append(x[0]);
	return matches;

def execute_command(command):
	# subprocess.call(command, shell=True);
	p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	[out, err] = p.communicate()
	returncode = p.returncode;
	return [out, err, returncode];

def measure_command(measure_file):
        return (SCRIPT_PATH + r'/../tools/cgmemtime/cgmemtime -o ' + measure_file + ' ');

def setup_measure_command():
	### Run a simple measure to try and see if the Cgroups was initialized.
	# command = '%s time' % (measure_command('/tmp/temp-measure.txt'));
	# [out, err, returncode] = execute_command(command);
	# if (returncode != 0 or ('Could not create new sub-cgroup' in out) or ('Could not create new sub-cgroup' in err)):
	# if (returncode != 0):
	sys.stderr.write('Sudo is needed for the Cgroups API (memory and time measurements).\n');
	#command = 'sudo %s/../tools/cgmemtime/cgmemtime --setup -g $USER --perm 775' % (SCRIPT_PATH);
	command = '%s/../tools/cgmemtime/cgmemtime --setup -g $USER --perm 775' % (SCRIPT_PATH);
	[out, err, returncode] = execute_command(command);
	# p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	# [out, err] = p.communicate()
	# if (('Could not create new sub-cgroup' in out) or ('Could not create new sub-cgroup' in err)):
	# 	sys.stderr.write('ERROR: Could not initialize Cgmemtime! Exiting.\n');
	# exit(1);



if __name__ == '__main__':
	pass;
