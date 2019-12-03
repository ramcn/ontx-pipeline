#! /usr/bin/python

import os;
import sys;
import utility_sam;

def parse_vcf_positions(vcf_file):
	try:
		fp = open(vcf_file, 'r');
		lines = fp.readlines();
		fp.close();
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % vcf_file);
		sys.stderr.write(str(e) + '\n');
		exit(1);

	positions = [];
	ref_bases = [];
	alt_bases = [];
	for line in lines:
		if (line[0] == '#'):
			continue;
		split_line = line.strip().split('\t');
		positions.append(int(split_line[1]) - 1);
		ref_bases.append(split_line[3]);
		alt_bases.append(split_line[4]);

	return [positions, ref_bases, alt_bases];



def verbose_usage_and_exit():
	sys.stderr.write('Takes a VCF file of simulated reference mutations, extracts positions, and checks alignments of simulated reads. For every position, number of correctly placed bases is counted.\n');
	sys.stderr.write('\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s <alignments.sam> <reference.sam> <mutations.vcf>' % sys.argv[0]);
	sys.stderr.write('\n');

	exit(0);

def main():
	if (len(sys.argv) != 4):
		verbose_usage_and_exit();

	query_sam = sys.argv[1];
	reference_sam = sys.argv[2];
	vcf_file = sys.argv[3];

	sys.stderr.write('Loading query SAM file...\n');
	[hashed_query, num_queries, num_unique_queries] = utility_sam.HashSAMWithFilter(query_sam, {});
	sys.stderr.write('Loading reference SAM file...\n');
	[hashed_reference, num_references, num_unique_references] = utility_sam.HashSAMWithFilter(reference_sam, {});
	sys.stderr.write('Loading positions from the VCF file...\n');
	[positions, ref_bases, alt_bases] = parse_vcf_positions(vcf_file);

	out_summary_prefix = os.path.splitext(vcf_file)[0];

	sys.stderr.write('Starting the counting process...\n');
	[accuracy, accuracy_called_bases] = utility_sam.CountCorrectlyMappedBasesAtPositions(hashed_query, hashed_reference, positions, ref_bases, alt_bases, out_summary_prefix=out_summary_prefix);
	sys.stderr.write('Accuracy: %.2f\n' % accuracy);
	sys.stderr.write('Accuracy (only called bases): %.2f\n' % accuracy_called_bases);



if __name__ == "__main__":
	main();