#! /usr/bin/python

SCRIPT_VERSION = 'v4';

import os;
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
import sys;

import re;
import math;
import time;
import numpy as np;

# import loadacc;
# from loadacc import AccuracyInfo;
from basicdefines import *;

# import main_sam_analysis;
import plot_with_seabourne_v2;
# import analyze_correctly_mapped_bases;



USE_MATPLOTLIB = True;
try:
	import matplotlib.pyplot as plt;
except:
	USE_MATPLOTLIB = False;

from matplotlib.font_manager import FontProperties;
import utility_sam;




def EvaluateAlignmentsFromPath(alignments_path, sam_suffix='', bp_dists=[10], count_correct_bases=True):
	dataset_name = alignments_path;
	current_folder_depth = len(alignments_path.split('/'));
	sam_files = find_files(alignments_path, '*%s.sam' % sam_suffix, (current_folder_depth));

	if len(sam_files) == 0:
		return None;

	split_out_path = alignments_path.split(EVALUATION_PATH_ROOT_ABS + '/');
	base_name = split_out_path[1];
	reference_sam = SCRIPT_PATH + '/../' + base_name + '/reads.sam';
	# out_scores_folder = SCRIPT_PATH + '/../' + RESULTS_ROOT + '/' + base_name;
	out_scores_folder = alignments_path + '';
	
	return EvaluateAlignments(reference_sam, sam_files, dataset_name, out_scores_folder, force_rerun=True, sam_suffix=sam_suffix, verbose_level=2, bp_dists=bp_dists, count_correct_bases=count_correct_bases);

### Example usage:
# EvaluateAlignments('reads-simulated/reads.sam', ['alignment/alignments.sam'], 'simulated_nanopore/escherichia_coli', 'alignments/');

def EvaluateAlignments(reference_sam, sam_files, dataset_name, out_scores_folder, force_rerun=False, sam_suffix='', verbose_level=0, plot_results=False, bp_dists=[10], count_correct_bases=True):
	filter_reads_with_multiple_alignments = False;
	reference_sequence_file = SCRIPT_PATH + '/../' + REFERENCE_GENOMES_ROOT + '/' + os.path.basename(dataset_name) + '.fa';
	
	if verbose_level > 0:
		print 'Dataset path:\n\t"%s"' % dataset_name;
		print 'Reference SAM:\n\t"%s"' % reference_sam;
		print 'SAM files:\n\t%s' % ('\n\t'.join(sam_files));
		print ' ';

	if not os.path.exists(out_scores_folder):
		print 'Creating output folders on path "%s".' % out_scores_folder;
		os.makedirs(out_scores_folder);
	
	qnames_with_multiple_alignments = {};
	if (filter_reads_with_multiple_alignments == True):
		print 'Finding multiple qname entries across all SAM files...';
		qnames_with_multiple_alignments = utility_sam.FindMultipleQnameEntries(sam_files);
		print 'Done finding multiple qname entries.';
		print ' ';
	else:
		print 'Using all reported alignments (even duplicates).';
	
	[hashed_reference, num_references, num_unique_references] = utility_sam.HashSAMWithFilter(reference_sam, qnames_with_multiple_alignments);

	out_folder_split = dataset_name.split('/');
	machine_name = out_folder_split[-2];
	genome_name = out_folder_split[-1];

#	total_prefix = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix;
	intermediate_results_path = out_scores_folder + '/analysis-intermediate'; # -' + machine_name + '-' + genome_name + '-' + sam_suffix;
	total_prefix_intermediate = intermediate_results_path + '/intermediate-' + sam_suffix;
	final_results_path = out_scores_folder + '/analysis-final'; # -' + machine_name + '-' + genome_name + '-' + sam_suffix;
	total_prefix = final_results_path + '/final-' + sam_suffix;

	if not os.path.exists(intermediate_results_path):
		print 'Creating output folders on path "%s".' % intermediate_results_path;
		os.makedirs(intermediate_results_path);
	if not os.path.exists(final_results_path):
		print 'Creating output folders on path "%s".' % final_results_path;
		os.makedirs(final_results_path);

	# print 'hashed_reference:';
	# print reference_sam;
	# print len(hashed_reference.keys());
	# for key in hashed_reference.keys():
	# 	# print hashed_reference[key][0].VerboseFormatLine();
	# 	# print hashed_reference[key][0].FormatAccuracy();
	# 	print hashed_reference[key][0].Verbose();
	# 	print '';
	# print '';

	ret_scores = {};

	# exit(1);

	all_accuracies = [];
	all_rmsd = [];
	# all_scores = [];
	all_out_path_prefixes = [];
	all_times = [];
	all_execution_times = [];
	
	all_distance_histograms = [];
	all_distance_histogram_percentage = [];
	all_mapq_histograms = [];
	all_mapq_histograms_percentage = [];
	all_labels = [];
	all_total_mapped = [];
	all_correctly_mapped_bases = [];
	all_correctly_mapped_bases_strict = [];
	all_correctly_mapped_bases_titles = [];
	
	total_sum_path = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix + '.sum';
	png_roc = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix + '-roc.png';
	png_precrec = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix + '-precision_recall.png';
	png_distances = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix + '-distance.png';
	png_mapq = out_scores_folder + '/total-' + machine_name + '-' + genome_name + '-' + sam_suffix + '-mapq.png';
	
	summary_line = '';
	all_summary_lines = [];
	all_roc_curves = [];
	all_precrec_curves = [];
	all_distance_accuracy = [];
	all_distance_recall = [];
	
	input_sam_path = '';
	if (len(sam_files) > 0):
		input_sam_path = os.path.basename(os.path.splitext(os.path.basename(sam_files[0]))[0]);
	
	sys.stderr.write('Counting the M operations in the input reference SAM:\n');
	# dataset_total_mapped_ref_num_m_ops = analyze_correctly_mapped_bases.CountMOps(hashed_reference);
	dataset_total_mapped_ref_num_m_ops = 3;

	current_sam_file = 0;
	for sam_file in sam_files:
		current_sam_file += 1;

		sys.stderr.write('Processing file: %s\n' % sam_file);
		sam_basename = os.path.splitext(os.path.basename(sam_file));

		if not os.path.exists(sam_file):
			sys.stderr.write('File "%s" does not exist! Skipping.\n' % sam_file);
			continue;
			# ret_scores[sam_basename[0]] = [distance_accuracy[1][10], distance_recall[1][10], percent_correctly_mapped_bases];

		out_scores_path_prefix = out_scores_folder + '/' + sam_basename[0];
		[is_sam_modified, modified_timestamp] = utility_sam.CheckSamModified(sam_file, out_scores_path_prefix);
		modified_time = time.ctime(float(modified_timestamp));
		
		execution_stats = utility_sam.GetExecutionStats(sam_file);
		
		sam_lines = ScoreAligned(hashed_reference, sam_file, shift_reverse_start=True);
		[hashed_sam_lines, num_sam_lines, num_unique_sam_lines] = utility_sam.HashSAMLines(sam_lines);
		[total_mapped, distance_histogram, distance_accuracy, distance_recall, mapq_histograms, mapq_histograms_percentage, summary_line] = CalculateStats(num_unique_references, sam_lines, sam_basename[0], execution_stats, modified_time, intermediate_results_path=intermediate_results_path);

		# print 'hashed_sam_lines:';
		# print reference_sam;
		# print len(hashed_sam_lines.keys());
		# for key in hashed_sam_lines.keys():
		# 	# print hashed_reference[key][0].VerboseFormatLine();
		# 	# print hashed_reference[key][0].FormatAccuracy();
		# 	print hashed_sam_lines[key][0].Verbose();
		# 	print '';
		# print '';

		all_distance_histograms.append(distance_histogram);
		all_distance_accuracy.append(distance_accuracy);
		all_distance_recall.append(distance_recall);
		all_mapq_histograms.append(mapq_histograms);
		all_mapq_histograms_percentage.append(mapq_histograms_percentage);
		all_total_mapped.append(total_mapped);
		all_labels.append(sam_basename);
		
		all_times.append(modified_time);
		all_out_path_prefixes.append(out_scores_path_prefix);
		all_execution_times.append(execution_stats);
		
		[roc, precrec] = GetROCFromEvaluatedSAM(sam_file, sam_lines, hashed_reference, 10, out_scores_path_prefix);
		all_roc_curves.append(roc);
		all_precrec_curves.append(precrec);

		# [percent_correctly_mapped_bases, num_correctly_mapped_bases, dataset_mapped_ref_num_m_ops] = analyze_correctly_mapped_bases.CountCorrectlyMappedBases(sam_file, hashed_reference, '');
		# [precision_correctly_mapped_bases, recall_correctly_mapped_bases, num_correctly_mapped_bases, dataset_mapped_ref_num_m_ops] = CountCorrectlyMappedBases(sam_basename[0], hashed_sam_lines, hashed_reference, '', use_strict=False);
		[precision_correctly_mapped_bases, recall_correctly_mapped_bases, num_correctly_mapped_bases, dataset_mapped_ref_num_m_ops] = [0.0, 0.0, 0, 0];
		[precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict, num_correctly_mapped_bases_strict, dataset_mapped_ref_num_m_ops_strict, num_ref_bases_strict] = [0.0, 0.0, 0, 0, 0];
		if (count_correct_bases == True):
			[precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict, num_correctly_mapped_bases_strict, dataset_mapped_ref_num_m_ops_strict, num_ref_bases_strict] = utility_sam.CountCorrectlyMappedBases(hashed_sam_lines, hashed_reference, '', sam_basename=sam_basename[0]);
		# exit(0);
			# [precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict, num_correctly_mapped_bases_strict, dataset_mapped_ref_num_m_ops_strict, num_ref_bases_strict] = utility_sam.CountCorrectlyMappedBasesAtPositions(hashed_sam_lines, hashed_reference, '', sam_basename=sam_basename[0]);
		# percent_correctly_mapped_bases = 0.0;
		# num_correctly_mapped_bases = 1;
		# dataset_mapped_ref_num_m_ops = 2;
# CountCorrectlyMappedBasesAtPositions

		# percent_correctly_mapped_bases_dataset_total = (float(num_correctly_mapped_bases) / float(dataset_total_mapped_ref_num_m_ops)) * 100.0;
		all_correctly_mapped_bases.append([precision_correctly_mapped_bases, recall_correctly_mapped_bases, num_correctly_mapped_bases, dataset_mapped_ref_num_m_ops, dataset_total_mapped_ref_num_m_ops]);
		all_correctly_mapped_bases_strict.append([precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict, num_correctly_mapped_bases_strict, dataset_mapped_ref_num_m_ops_strict, dataset_total_mapped_ref_num_m_ops]);
		all_correctly_mapped_bases_titles = ['precision_correctly_mapped_bases', 'recall_correctly_mapped_bases', 'num_correctly_mapped_bases', 'dataset_mapped_ref_num_m_ops', 'dataset_total_mapped_ref_num_m_ops'];
		summary_line_correct_bases = '';
		# summary_line_correct_bases += 'Percent correctly mapped bases: %.2f\n' % percent_correctly_mapped_bases;
		# summary_line_correct_bases += 'Number of correctly mapped bases: %d\n' % num_correctly_mapped_bases;
		# summary_line_correct_bases += 'Number of M ops in mapped reads: %d\n' % dataset_mapped_ref_num_m_ops;
		# summary_line_correct_bases += 'Number of M ops in all input reads: %d\n' % dataset_total_mapped_ref_num_m_ops;
		summary_line_correct_bases += 'Number of correctly mapped bases: %d\n' % num_correctly_mapped_bases_strict;
		summary_line_correct_bases += 'Number of mapped bases: %d\n' % dataset_mapped_ref_num_m_ops_strict;
		summary_line_correct_bases += 'Number of bases in input SAM: %d\n' % num_ref_bases_strict;

		summary_line_correct_bases += 'Precision (per-base): %.2f\n' % (precision_correctly_mapped_bases);
		summary_line_correct_bases += 'Recall (per-base): %.2f\n' % (recall_correctly_mapped_bases);
		summary_line_correct_bases += '(strict) Precision (per-base): %.2f\n' % (precision_correctly_mapped_bases_strict);
		summary_line_correct_bases += '(strict) Recall (per-base): %.2f\n' % (recall_correctly_mapped_bases_strict);
		# summary_line_correct_bases += 'Percent of correct M ops in all input reads: %.2f\n' % ((float(num_correctly_mapped_bases) / float(dataset_total_mapped_ref_num_m_ops)) * 100.0);
		summary_line += summary_line_correct_bases;

		all_summary_lines.append(summary_line);

		

		sys.stdout.write(summary_line + '\n');


		# ret_scores[sam_basename[0]] = [	distance_accuracy[1][bp_dist], distance_recall[1][bp_dist],
		# 								precision_correctly_mapped_bases, recall_correctly_mapped_bases,
		# 								precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict];
		ret_scores[sam_basename[0]] = [];
		for bp_dist in bp_dists:
			bp_precrec = [(bp_dist, distance_accuracy[1][bp_dist], distance_recall[1][bp_dist])];
			ret_scores[sam_basename[0]] += bp_precrec;

		ret_scores[sam_basename[0]] += [(0, precision_correctly_mapped_bases, recall_correctly_mapped_bases),
										(0, precision_correctly_mapped_bases_strict, recall_correctly_mapped_bases_strict)];

		[cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit] = utility_sam.ParseMemTime(sam_file);
		ret_scores[sam_basename[0]] += [(time_unit, 0, cputime)];
		ret_scores[sam_basename[0]] += [(mem_unit, 0, maxrss)];


		
		accuracies_path = intermediate_results_path + '/' + sam_basename[0];
		if verbose_level > 0:
			print 'Writing accuracies to prefix: "%s"' % total_prefix_intermediate;
		WriteAccuracies(hashed_sam_lines, accuracies_path, write_duplicates=False);

		WriteSummary(sam_files[0:current_sam_file], all_summary_lines, total_prefix + '.sum');
		WriteHistogramsDistance(input_sam_path, all_distance_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance.hist' % (total_prefix));
		WriteHistogramsDistance(input_sam_path, all_distance_accuracy, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance_correct_mapped.csv' % (total_prefix));
		WriteHistogramsDistance(input_sam_path, all_distance_recall, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance_correct_w_unmapped.csv' % (total_prefix));
		WriteHistogramsMapq(input_sam_path, all_mapq_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-mapq.hist' % (total_prefix));
		
		WriteROC(input_sam_path, all_roc_curves, all_labels, machine_name, genome_name, '%s-roc.csv' % (total_prefix));
		WriteROC(input_sam_path, all_precrec_curves, all_labels, machine_name, genome_name, '%s-precision_recall.csv' % (total_prefix));
	
		WriteCorrectBaseCount(input_sam_path, all_correctly_mapped_bases_titles, all_correctly_mapped_bases, all_labels, machine_name, genome_name , '%s-correct_bases.csv' % (total_prefix));
		WriteCorrectBaseCount(input_sam_path, all_correctly_mapped_bases_titles, all_correctly_mapped_bases_strict, all_labels, machine_name, genome_name , '%s-strict-correct_bases.csv' % (total_prefix));

	if verbose_level > 0:
		print ' ';
		print '====================================';
	
	### PlotHistogramsDistance(all_distance_recall, all_labels, machine_name + ', ' + genome_name, png_distances);
	#PlotHistogramsMapq(all_mapq_histograms_percentage, all_labels);
	#PlotROC(all_roc_curves, all_labels, 'Fraction of reads [%]', 'Correctly mapped reads [%]', 'Curve of mapping accuracy (distance 0)', machine_name + ', ' + genome_name, png_roc);
	### PlotROC(all_roc_curves, all_labels, 'False positive rate', 'True positive rate', 'Curve of mapping accuracy (distance 10)', machine_name + ', ' + genome_name, png_roc);
	### PlotROC(all_precrec_curves, all_labels, 'Recall', 'Precision', 'Precision-Recall wrt. MAPQ or AS (allowed distance 10)', machine_name + ', ' + genome_name, png_precrec);
	
	# WriteSummary(sam_files, all_summary_lines, total_prefix + '.sum');
	# WriteHistogramsDistance(input_sam_path, all_distance_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance.hist' % (total_prefix));
	# WriteHistogramsDistance(input_sam_path, all_distance_accuracy, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance_correct_mapped.csv' % (total_prefix));
	# WriteHistogramsDistance(input_sam_path, all_distance_recall, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-distance_correct_w_unmapped.csv' % (total_prefix));
	# WriteHistogramsMapq(input_sam_path, all_mapq_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, '%s-mapq.hist' % (total_prefix));
	
	# WriteROC(input_sam_path, all_roc_curves, all_labels, machine_name, genome_name, '%s-roc.csv' % (total_prefix));
	# WriteROC(input_sam_path, all_precrec_curves, all_labels, machine_name, genome_name, '%s-precision_recall.csv' % (total_prefix));
	
	if (plot_results == True):
		if verbose_level > 0:
			print 'Plotting simulated results...';
		plot_with_seabourne_v2.PlotSimulatedResults(total_prefix, sam_suffix, genome_name, SCRIPT_VERSION);

	if verbose_level > 0:
		print 'Post-analyzing SAM files...';
	
	#main_sam_analysis.ProcessSamFiles(sam_files, '-', out_scores_folder, num_unique_references, suppress_error_messages=True);
#	main_sam_analysis.ProcessSamFiles(sam_files, reference_sequence_file, out_scores_folder, num_unique_references, suppress_error_messages=True);
	
	if verbose_level > 0:
		print '====================================';
	
	#plt.show();
	
	return ret_scores;

def IncreaseHashCounts(hash_counts, param_name, param_increase_count):
	try:
		hash_counts[param_name] += param_increase_count;
	except:
		hash_counts[param_name] = param_increase_count;

# Scores the alignments in a single SAM filed, compared to the ground-truth. The ground-truth is given through the
# 'hashed_reference' parameter, obtained with the HashSAM function. The 'num_references' parameter is the total number of
# reads that were in the input dataset (can also be obtained from the HashSAM function). It is used for determining the
# limits for the x and y axis of ROC curves (although not used in this function, it is only passed through to the return
# value).
def ScoreAligned(hashed_reference, aligned_path, shift_reverse_start=True):
	sam_lines = [];
	
	try:
		fp_aligned = open(aligned_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fp_aligned;
		return [correct, wrong, num_references, roc_fp, roc_tp];
	
	for value in hashed_reference.values():
		for reference in value:
			reference.evaluated = 0;

	sam_basename = os.path.splitext(os.path.basename(aligned_path))[0];
	
  
 
	numDeformedHeader = 0;
	num_duplicates = 0;
	qname_counts = {};
	
	num_lines = 0;
	for line in fp_aligned:
		line = line.strip();
		
		if len(line) == 0:
			continue;
      
		if line[0] == '@':
			continue;
		
		is_orientation_and_chrom_ok = 0;
		hit_reference_pos = -1;
		hit_mapped_pos = -1;
		hit_min_distance = -1;
		hit_ref_reverse = False;
		is_header_deformed = False
		
		sam_line = utility_sam.SAMLine(line, sam_basename);
		header_key = utility_sam.TrimHeader(sam_line.qname);
		
		IncreaseHashCounts(qname_counts, header_key, 1);
		
		ref_header = '*';
		map_ref_header = sam_line.rname;
		
		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
		# uniformly evaluate the results!
		if 'blasr' in aligned_path.lower():
			#if sam_line.IsReverse():
				#sam_line.pos += sam_line.isize;
				#sam_line.clipped_pos += sam_line.isize;
			header_key = '/'.join(header_key.split('/')[:-1]);
			if sam_line.clip_count_front != 0 or sam_line.clip_count_back != 0:
				print 'BLASR CIGAR contains clipping! Please revise clipped_pos! Read: "%s".' % sam_line.qname;
		
		sam_line.trimmed_header = header_key;
		
		try:
			#m130404_014004_sidney_c100506902550000001823076808221337_s1_p0/25/0_593
#			href = hashed_reference[header_key.split('/')[0]];			# TODO: Ovo sam zadnje promijenio, radilo je dobro sa simuliranim podatcima, ali s pravim readovima je krsilo sve, jer su u headerima imali slasheve! Promjena napravljena 27.08.2014.
			href = hashed_reference[header_key];
			is_href_duplicate = False;
			for reference in href:
				if (reference.is_filtered_out == True):
					is_href_duplicate = True;
					break;
			
			if (is_href_duplicate == False):
				CompareAlignments(href, sam_line, shift_reverse_start);
			else:
				#num_duplicates += sam_line.is_duplicate;
				num_duplicates += 1;
				sam_line.is_duplicate = 1;
			
		except Exception, e:
			numDeformedHeader += 1;
			sam_line.is_correct_ref_and_orient = 0;
			is_header_deformed = True;
			sam_line.is_header_deformed = 1;
			pass;

		sam_lines.append(sam_line);

		num_lines += 1;

	fp_aligned.close();
	
	for sam_line in sam_lines:
		header_key = utility_sam.TrimHeader(sam_line.qname);
		sam_line.num_occurances_in_sam_file = qname_counts[header_key];
	
	print 'numDeformedHeader = %d' % numDeformedHeader;
	print 'numDuplicates = %d' % num_duplicates;
	
	return sam_lines;

# Compares a single SAM entry to a list of 'correct' alignments. Here, 'correct' alignments are given as a list,
# as forward and reverse SAM entries can have the same QNAME fields. Also, this allows for comparison with alternative
# alignments, for instance in case of real data that can be mapped to several different locations on the reference genome.
def CompareAlignments(reference_alignments, query, shift_reverse_start=True):
	FLAG_REVERSE = 1 << 4;			# Sequence is reversed.
	
	ret_reference_pos = -1;
	ret_mapped_pos = -1;
	ret_min_distance = -1;
	ret_hit = 0;
	ret_ref_reverse = False;
	ret_ref_header = '*';
	
	query.is_correct_ref_and_orient = 0;
	query.is_header_deformed = False;
	
	i = 0;
	for reference in reference_alignments:
		query_ref_name = query.rname;
		reference_ref_name = reference.rname;
		
		query_pos = query.clipped_pos;
		reference_pos = reference.clipped_pos;
		
		query_orientation = query.IsReverse();
		reference_orientation = reference.IsReverse();

		# Distance from the mapped position to the current possible location (multiple
		# locations are taken into account).
		distance = abs(reference_pos - query_pos);
		
		# Track the minimum distance and the exact positions for those. This is needed
		# even if the read is not correctly mapped, because we want to see where it
		# ended up.
		if ret_min_distance == -1 or distance < ret_min_distance:
			ret_min_distance = distance;
			ret_reference_pos = reference_pos;
			ret_mapped_pos = query_pos;
			ret_ref_reverse = reference.IsReverse();
			ret_ref_header = reference_ref_name;

		# Here we check if the alignment's orientation  is valid.
		if ((query_ref_name == reference_ref_name or
		     query_ref_name.startswith(reference_ref_name) or
		     reference_ref_name.startswith(query_ref_name)) and
		    (query_orientation == reference_orientation)):
			if reference.evaluated > 0:
				#print '\rRead "%s" has already been evaluated! Is this an alternative position?' % query.qname;
				if (reference.min_distance >= distance):
					reference.min_distance = distance;
				query.is_duplicate = 1;
			else:
				reference.min_distance = distance;
				
			reference.evaluated += 1;
			ret_hit = 1;
			query.is_correct_ref_and_orient = 1;
			
		i += 1;

	query.min_distance = ret_min_distance;
	query.actual_ref_pos = ret_reference_pos;
	query.actual_ref_reverse = ret_ref_reverse;
	query.actual_ref_header = ret_ref_header;
	query.mapped_pos_with_shift = ret_mapped_pos;

def GetROCFromEvaluatedSAM(sam_file, sam_lines, hashed_reference, allowed_distance, out_path_prefix):
	#sorted_lines_by_quality = sorted(sam_lines, reverse=True, key=lambda sam_line: (sam_line.chosen_quality, -sam_line.min_distance));
	sorted_lines_by_quality = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.chosen_quality);
	
	#fp = open('temp.txt', 'w');
	#fp.write('\n'.join([sam_line.FormatAccuracy() for sam_line in sorted_lines_by_quality]));
	#fp.close();
	
	total_tp = 0.0;
	total_fp = 0.0;
	total_fn = 0.0;
	total_tn = 0.0;
	total_not_mapped = 0.0;
	tp_list = [];
	fp_list = [];
	fn_list = [];
	tn_list = [];
	#total_not_mapped_list = [];
	
	tpr_list = [];
	fpr_list = [];
	
	precision_list = [];
	recall_list = [];
	
	num_tp = 0;
	num_fp = 0;
	
	unique_reads = {};
	
	max_mapq = 0;
	
	current_num_alignments = 0;
	num_mapped_alignments = 0;
	for sam_line in sorted_lines_by_quality:
		current_num_alignments += 1;
		unique_reads[sam_line.qname] = 1;
		
		if (sam_line.IsMapped() == True):
			num_mapped_alignments += 1;
			if (sam_line.chosen_quality > max_mapq):
				max_mapq = sam_line.chosen_quality + 0.0;
		
			if (sam_line.min_distance <= allowed_distance):
				total_tp += (1.0 / sam_line.num_occurances_in_sam_file);
				num_tp += 1;
			else:
				total_fp += (1.0 / sam_line.num_occurances_in_sam_file);
				num_fp += 1;
		else:
			total_not_mapped += (1.0 / sam_line.num_occurances_in_sam_file);
	
	total_not_mapped += (len(hashed_reference.keys()) - len(unique_reads.keys()));
	
	tp = 0.0;
	fp = 0.0;
	fn = 0.0;
	tn = 0.0;
	
	i = 0;
	i_mapped = 0;
	previous_mapq = 0;
	while i < len(sorted_lines_by_quality):
		#print 'i1 = %d' % i;
		sam_line = sorted_lines_by_quality[i];
		
		if (sam_line.IsMapped() == False):
			i += 1;
			continue;
		
		i_mapped += 1;
		
		current_mapq = sam_line.chosen_quality + 0.0;
		if (i == 0):
			previous_mapq = current_mapq + 0.0;
		
		#print 'i2 = %d' % i;
		
		if (current_mapq != previous_mapq and (i > 0 and sorted_lines_by_quality[i-1].IsMapped()==True)):
			fn = total_tp - tp + total_not_mapped;
			tn = total_fp - fp + total_not_mapped;
	
			# print 'tp = %d' % tp;
			# print 'fp = %d' % fp;
			# print 'fn = %d' % fn;

			precision = float(tp) / (tp + fp) if ((tp + fp) != 0) else 0;
			recall = float(tp) / (tp + fn) if ((tp + fn) != 0) else 0;
			precision_list.append(precision);
			recall_list.append(recall);
			#tpr_list.append(float(tp) / float(tp + fn));
			#fpr_list.append(float(fp) / float(fp + tn));
			#fpr_list.append(float(i) / float(num_mapped_alignments));
			#tpr_list.append(float(tp) / float(total_tp + total_fp));
			#fpr_list.append(1.0 - float(previous_mapq) / float(max_mapq));
			tpr_list.append(float(tp) / float(current_num_alignments));
			fpr_list.append(float(fp) / float(current_num_alignments));
			
		if (sam_line.min_distance <= allowed_distance):
			tp += (1.0 / sam_line.num_occurances_in_sam_file);
		else:
			fp += (1.0 / sam_line.num_occurances_in_sam_file);
		
		previous_mapq = current_mapq + 0.0;
		
		i += 1;
	
	fn = total_tp - tp + total_not_mapped;

	precision = tp / (tp + fp) if ((tp + fp != 0)) else 0.0;
	recall = tp / (tp + fn) if ((tp + fn) != 0) else 0.0;
	precision_list.append(precision);
	recall_list.append(recall);
	tpr_list.append(float(tp) / float(current_num_alignments) if (current_num_alignments != 0) else 0.0);
	fpr_list.append(float(fp) / float(current_num_alignments) if (current_num_alignments != 0) else 0.0);

	roc = [fpr_list, tpr_list];
	precrec = [recall_list, precision_list];
	
	return [roc, precrec];

	#tpr_list.append(tp / total_tp);
	#fpr_list.append(fp / total_fp);
	#fpr_list.append(float(i) / float(num_mapped_alignments));
	#tpr_list.append(float(tp) / float(total_tp + total_fp));
	#fpr_list.append(1.0 - float(previous_mapq) / float(max_mapq));
	
	#tpr_list = [(value / float(tp)) for value in tp_list];
	#fpr_list = [(value / float(fp)) for value in fp_list];
	##precision_list = [(value / float(num_mapped_alignments)) for value in tp_list];
	##recall_list = [(value / float(len(unique_reads.keys()))) for value in tp_list];
	#precision_list = [(value / float(tp + fp)) for value in tp_list];
	#recall_list = [(value / float(tp + fn)) for value in tp_list];
	#print precision_list;
	
	
	
	#Precision = TP / (TP+FP)
	#Recall = TP / (TP+FN)
	#Accuracy = (TP + TN) / (TP + TN + FP + FN)

	#tpr_list = [(float(tp) / num_tp) for tp in tp_list];
	#fpr_list = [(float(fp) / num_fp) for fp in fp_list];
	
	#roc_curves = [[tpr_list, fpr_list]];
	#labels = [os.path.splitext(os.path.basename(sam_file))[0]];
	#dataset_description = sam_file;
	#out_png_path = '';
	#PlotROC(roc_curves, labels, dataset_description, out_png_path);
	#plt.show();
	
	#fp_out = open(out_path_prefix + '_v3.roc', 'w');
	##i = 0;
	##while i < len(tpr_list):
		##fp_out.write('%d\t%d\n' % (tpr_list[i], fpr_list[i]));
		##i += 1;
	#i = 0;
	#while i < len(sorted_lines_by_quality):
		#fp_out.write('%s\n' % (sorted_lines_by_quality[i].FormatAccuracy()));
		#i += 1;
	#fp_out.close();
	
	#print 

	#return [tp_list, fp_list, tpr_list, fpr_list];
	#return [tpr_list, fpr_list];

def CountSameCigarOps(tested_cigar_pos_list, reference_cigar_pos_list):
	num_same_ops = 0;
	num_same_m_ops = 0;
	num_same_i_ops = 0;
	num_same_d_ops = 0;
	num_same_other_ops = 0;
	
	i = 0;
	j = 0;
	
	while (i < len(tested_cigar_pos_list) and j < len(reference_cigar_pos_list)):
		# If the coordinates are not equal, increase one of them and continue.
		if (tested_cigar_pos_list[i][2] < reference_cigar_pos_list[j][2]):
			i += 1;
			continue;
		elif (tested_cigar_pos_list[i][2] > reference_cigar_pos_list[j][2]):
			j += 1;
			continue;
		# If the coordinates are equal, then compare the CIGAR operations.
		else:
			# If the operation and the count are the same, we found a match.
			if (tested_cigar_pos_list[i][0] == reference_cigar_pos_list[j][0] and
			    tested_cigar_pos_list[i][1] == reference_cigar_pos_list[j][1]):
				num_same_ops += 1;
				
				# This counts only match/mismatch operations that are the same.
				if (tested_cigar_pos_list[i][1] in 'M=X'):
					num_same_m_ops += 1;
				# This counts only insertion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'I'):
					num_same_i_ops += 1;
				# This counts only deletion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'D'):
					num_same_d_ops += 1;
				# This counts all other operations that are the same.
				else:
					num_same_other_ops += 1;
					
				#print tested_cigar_pos_list[i];
			else:
				if (tested_cigar_pos_list[i][1] == 'I' and reference_cigar_pos_list[j][1] != 'I'):
					i += 1;
					continue;
				elif (tested_cigar_pos_list[i][1] != 'I' and reference_cigar_pos_list[j][1] =='I'):
					j += 1;
					continue;
				
				#i = (i + 1) if (tested_cigar_pos_list[i][1] == 'I') else i;
				#j = (j + 1) if (reference_cigar_pos_list[j][1] == 'I') else j;
				#continue;
			
		i += 1;
		j += 1;
		
	return [num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops];

def CountSameCigarOpsStrict(tested_cigar_pos_list, reference_cigar_pos_list):
	num_same_ops = 0;
	num_same_m_ops = 0;
	num_same_i_ops = 0;
	num_same_d_ops = 0;
	num_same_other_ops = 0;
	
	i = 0;
	j = 0;
	
	while (i < len(tested_cigar_pos_list) and j < len(reference_cigar_pos_list)):
		# If the coordinates are not equal, increase one of them and continue.
		if (tested_cigar_pos_list[i][2] < reference_cigar_pos_list[j][2]):
			i += 1;
			continue;
		elif (tested_cigar_pos_list[i][2] > reference_cigar_pos_list[j][2]):
			j += 1;
			continue;
		# If the coordinates are equal, then compare the CIGAR operations.
		else:
			# If the operation and the count are the same, we found a match.
			if (tested_cigar_pos_list[i][0] == reference_cigar_pos_list[j][0] and
			    tested_cigar_pos_list[i][1] == reference_cigar_pos_list[j][1] and
			    tested_cigar_pos_list[i][3] == reference_cigar_pos_list[j][3]):
				num_same_ops += 1;
				
				# This counts only match/mismatch operations that are the same.
				if (tested_cigar_pos_list[i][1] in 'M=X'):
					num_same_m_ops += 1;
				# This counts only insertion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'I'):
					num_same_i_ops += 1;
				# This counts only deletion operations that are the same.
				elif (tested_cigar_pos_list[i][1] == 'D'):
					num_same_d_ops += 1;
				# This counts all other operations that are the same.
				else:
					num_same_other_ops += 1;
					
				#print tested_cigar_pos_list[i];
			else:
				if (tested_cigar_pos_list[i][1] == 'I' and reference_cigar_pos_list[j][1] != 'I'):
					i += 1;
					continue;
				elif (tested_cigar_pos_list[i][1] != 'I' and reference_cigar_pos_list[j][1] =='I'):
					j += 1;
					continue;
				
				#i = (i + 1) if (tested_cigar_pos_list[i][1] == 'I') else i;
				#j = (j + 1) if (reference_cigar_pos_list[j][1] == 'I') else j;
				#continue;
			
		i += 1;
		j += 1;
		
	return [num_same_ops, num_same_m_ops, num_same_i_ops, num_same_d_ops, num_same_other_ops];

def CountOperations(cigar_list):
	num_m_ops = 0;
	num_i_ops = 0;
	num_d_ops = 0;
	num_other_ops = 0;
	
	for cigar in cigar_list:
		if (cigar[1] in 'M=X'):
			num_m_ops += 1;
		elif (cigar[1] == 'I'):
			num_i_ops += 1;
		elif (cigar[1] == 'D'):
			num_d_ops += 1;
		else:
			num_other_ops += 1;

	num_m_ops = 1 if (num_m_ops == 0) else num_m_ops;
	num_i_ops = 1 if (num_i_ops == 0) else num_i_ops;
	num_d_ops = 1 if (num_d_ops == 0) else num_d_ops;
	num_other_ops = 1 if (num_other_ops == 0) else num_other_ops;
	
	return [num_m_ops, num_i_ops, num_d_ops, num_other_ops];

# def CompareBasePositions(query_sam, ref_sam):
# 	qsam_ref_coords = [None] * query_sam.CalcReadLengthFromCigar();
# 	rsam_ref_coords = [None] * ref_sam.CalcReadLengthFromCigar();

# 	num_mapped_bases = len(qsam_ref_coords) - query_sam.clip_count_front - query_sam.clip_count_back;
# 	num_ref_bases = len(rsam_ref_coords) - ref_sam.clip_count_front - ref_sam.clip_count_back;

# 	if (len(qsam_ref_coords) < len(rsam_ref_coords)):
# 		# sys.stderr.write('Warning: Mappers output does not conform to SAM format specification! CIGAR field does not specify sequence of equal length as in the input FASTA file. Possibly hard clipping operations are missing.\n');	
# 		return [0, num_mapped_bases, num_ref_bases];

# 	query_cigpos = query_sam.CalcCigarStartingPositions(False);
# 	for cigpos in query_cigpos:
# 		[cig_count, cig_op, pos_on_ref, pos_on_query] = cigpos;
# 		if (cig_op in 'M=X'):
# 			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = range(pos_on_ref, (pos_on_ref + cig_count));
# 		elif (cig_op == 'I'):
# 			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-1]*cig_count;
# 		elif (cig_op == 'S'):
# 			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-3]*cig_count;
# 		elif (cig_op == 'H'):
# 			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-4]*cig_count;

# 	ref_cigpos = ref_sam.CalcCigarStartingPositions(False);
# 	for cigpos in ref_cigpos:
# 		[cig_count, cig_op, pos_on_ref, pos_on_query] = cigpos;
# 		if (cig_op in 'M=X'):
# 			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = range(pos_on_ref, (pos_on_ref + cig_count));
# 		elif (cig_op == 'I'):
# 			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-1]*cig_count;
# 		elif (cig_op == 'S'):
# 			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-3]*cig_count;
# 		elif (cig_op == 'H'):
# 			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-4]*cig_count;

# 	num_correct_bases = 0;
# 	i = 0;
# 	while (i < len(qsam_ref_coords)):
# 		# Skip the clipped bases:
# 		if (qsam_ref_coords[i] != -3 and qsam_ref_coords[i]!= -4):
# 			# Count the equal bases
# 			if (qsam_ref_coords[i] == rsam_ref_coords[i]):
# 				num_correct_bases += 1;
# 			# else:
# 			# 	# Just debug output.
# 			# 	print 'qsam[i] = %d\trsam[i] = %d\tseq[i] = %s' % (qsam_ref_coords[i], rsam_ref_coords[i], ref_sam.seq[i]);
# 		i += 1;

# 	# print query_sam.cigar;
# 	# print qsam_ref_coords;
# 	# print query_sam.original_line;
# 	# print '';
# 	# print 'num_correct_bases = %d' % num_correct_bases;
# 	# print 'num_mapped_bases = %d' % num_mapped_bases;
# 	# print ref_sam.original_line;

# 	return [num_correct_bases, num_mapped_bases, num_ref_bases];



	# exit(1);

# def CountCorrectlyMappedBases(sam_basename, hashed_sam_lines, hashed_reference_sam, out_summary_prefix='', use_strict=False):
# 	if (use_strict == False):
# 		return [0.0, 0.0, 0, 0];

# 	fp_out = None;
# 	out_file = out_summary_prefix + '.csv';
# 	if (out_summary_prefix != ''):
# 		try:
# 			fp_out = open(out_file, 'w');
# 		except IOError:
# 			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
# 			exit(1);
# 	sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');


# 	total_ref_bases = 0;
# 	for qname in hashed_reference_sam.keys():
# 		ref_sam = hashed_reference_sam[qname][0];
# 		# total_ref_bases += (ref_sam.CalcReadLengthFromCigar() - ref_sam.clip_count_front - ref_sam.clip_count_back);
# 		num_ref_bases = len(ref_sam.seq);
# 		# Remove the counts of soft-clipped bases from the read length. We ignore hard clipped bases, because they are
# 		# not present in the SEQ field anyway.
# 		num_ref_bases -= (ref_sam.clip_count_front if ref_sam.clip_op_front == 'S' else 0);
# 		num_ref_bases -= (ref_sam.clip_count_back if ref_sam.clip_op_back == 'S' else 0);
# 		total_ref_bases += num_ref_bases;

# 	sum_correct_bases = 0;
# 	sum_mapped_bases = 0;
# 	sum_ref_bases = 0;

# 	i = 0;
# 	for qname in hashed_sam_lines.keys():
# 		i += 1;
# 		if ((i % 100) == 0):
# 			sys.stderr.write('\rLine %d' % (i));
# 			sys.stderr.flush();

# 		sam_line = hashed_sam_lines[qname][0];

# 		if (sam_line.IsMapped() == False):
# 			continue;
# 		if ((qname in hashed_sam_lines) == False):
# 			sys.stderr.write('\tERROR: Reference SAM does not contain qname "%s"!\n' % (qname));
# 			continue;

# 		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
# 		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
# 		# uniformly evaluate the results!
# 		if 'blasr' in sam_basename.lower():
# 			qname = '/'.join(qname.split('/')[:-1]);
# 			if sam_line.clip_count_front != 0 or sam_line.clip_count_back != 0:
# 				print 'BLASR CIGAR contains clipping! Please revise clipped_pos! Read: "%s".' % sam_line.qname;

# 		sam_reference = hashed_reference_sam[qname][0];

# 		[num_correct_bases, num_mapped_bases, num_ref_bases] = CompareBasePositions(sam_line, sam_reference);
# 		sum_correct_bases += num_correct_bases;
# 		sum_mapped_bases += num_mapped_bases;
# 		sum_ref_bases += num_ref_bases;

# 		if (float(num_correct_bases) / float(num_mapped_bases) < 0.75):
# 			print 'Original line:';
# 			print sam_line.original_line;
# 			print '';
# 			print 'Reference SAM line:';
# 			print sam_reference.original_line;
# 			print '';
# 			print '';

# 		# Trebao bih ovdje dohvatiti broj tocno mapiranih baza, i broj baza koje ukupno postoje, i onda izracunati statistike i vratiti ih dolje!

# 	precision = 100.0 * float(sum_correct_bases) / float(sum_mapped_bases);
# 	recall = 100.0 * float(sum_correct_bases) / float(total_ref_bases);

# 	if (out_summary_prefix != ''):
# 		fp_out.write('percent_correct_m\tnum_correct_m\tnum_m_ops_in_reference\n');
# 		fp_out.write('%.2f\t%.2f\t%.2f\t%.2f\n' % (precision, recall, sum_correct_bases, sum_mapped_bases));
# 		fp_out.close();
# 	sys.stderr.write('\n');

# 	return [precision, recall, sum_correct_bases, sum_mapped_bases, total_ref_bases];



# def CountCorrectlyMappedBases2(sam_basename, hashed_sam_lines, hashed_reference_sam, out_summary_prefix='', use_strict=False):
# 	fp_out = None;
	
# 	out_file = out_summary_prefix + '.csv';
	
# 	if (out_summary_prefix != ''):
# 		try:
# 			fp_out = open(out_file, 'w');
# 		except IOError:
# 			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
# 			exit(1);

# 	sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');

# 	total_num_same_ops_individual_m = 0;
# 	total_num_same_m_ops_individual_m = 0;
# 	total_num_same_i_ops_individual_m = 0;
# 	total_num_same_d_ops_individual_m = 0;
# 	total_num_same_other_ops_individual_m = 0;
# 	query_num_m_ops_mapped_reads = 0;
# 	ref_num_m_ops_mapped_reads = 0;
# 	ref_num_i_ops_mapped_reads = 0;
# 	ref_num_d_ops_mapped_reads = 0;
# 	dataset_total_ref_num_other_ops = 0;

# 	ref_num_m_ops_all_reads = 0;
# 	ref_num_i_ops_all_reads = 0;
# 	ref_num_d_ops_all_reads = 0;

# 	reference_cigar_pos_list_individual_m = {};
# 	for qname in hashed_reference_sam.keys():
# 		sam_reference = hashed_reference_sam[qname][0];
# 		cigar_pos_list = sam_reference.CalcCigarStartingPositions(True);
# 		reference_cigar_pos_list_individual_m[qname] = cigar_pos_list;
# 		[ref_num_m_ops_individual_m, ref_num_i_ops_individual_m, ref_num_d_ops_individual_m, ref_num_other_ops_individual_m] = CountOperations(cigar_pos_list);
# 		ref_num_m_ops_all_reads += ref_num_m_ops_individual_m;
# 		ref_num_i_ops_all_reads += ref_num_i_ops_individual_m;
# 		ref_num_d_ops_all_reads += ref_num_d_ops_individual_m;

# 	i = 0;
# 	for qname in hashed_sam_lines.keys():
# 		i += 1;
# 		if ((i % 100) == 0):
# 			sys.stderr.write('\rLine %d' % (i));
# 			sys.stderr.flush();

# 		sam_line = hashed_sam_lines[qname][0];

# 		if (sam_line.IsMapped() == False):
# 			continue;
# 		if ((qname in hashed_sam_lines) == False):
# 			sys.stderr.write('\tERROR: Reference SAM does not contain qname "%s"!\n' % (qname));
# 			continue;

# 		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
# 		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
# 		# uniformly evaluate the results!
# 		if 'blasr' in sam_basename.lower():
# 			qname = '/'.join(qname.split('/')[:-1]);
# 			if sam_line.clip_count_front != 0 or sam_line.clip_count_back != 0:
# 				print 'BLASR CIGAR contains clipping! Please revise clipped_pos! Read: "%s".' % sam_line.qname;

# 		sam_reference = hashed_reference_sam[qname][0];

# 		cigar_pos_list_individual_m = sam_line.CalcCigarStartingPositions(True);
# 		# reference_cigar_pos_list_individual_m = sam_reference.CalcCigarStartingPositions(True);
# 		reference_cigar_pos_list = reference_cigar_pos_list_individual_m[qname];
		
# 		[num_m_ops_individual_m, num_i_ops_individual_m, num_d_ops_individual_m, num_other_ops_individual_m] = CountOperations(cigar_pos_list_individual_m);
# 		[ref_num_m_ops_individual_m, ref_num_i_ops_individual_m, ref_num_d_ops_individual_m, ref_num_other_ops_individual_m] = CountOperations(reference_cigar_pos_list);
		
# 		query_num_m_ops_mapped_reads += num_m_ops_individual_m;
# 		ref_num_m_ops_mapped_reads += ref_num_m_ops_individual_m;
# 		ref_num_i_ops_mapped_reads += ref_num_i_ops_individual_m;
# 		ref_num_d_ops_mapped_reads += ref_num_d_ops_individual_m;
# 		dataset_total_ref_num_other_ops += ref_num_other_ops_individual_m;

# 		if (use_strict == False):
# 			[num_same_ops_individual_m, num_same_m_ops_individual_m, num_same_i_ops_individual_m, num_same_d_ops_individual_m, num_same_other_ops_individual_m] = CountSameCigarOps(cigar_pos_list_individual_m, reference_cigar_pos_list);
# 		else:
# 			[num_same_ops_individual_m, num_same_m_ops_individual_m, num_same_i_ops_individual_m, num_same_d_ops_individual_m, num_same_other_ops_individual_m] = CountSameCigarOpsStrict(cigar_pos_list_individual_m, reference_cigar_pos_list);
# 		total_num_same_ops_individual_m += num_same_ops_individual_m;
# 		total_num_same_m_ops_individual_m += num_same_m_ops_individual_m;
# 		total_num_same_i_ops_individual_m += num_same_i_ops_individual_m;
# 		total_num_same_d_ops_individual_m += num_same_d_ops_individual_m;
# 		total_num_same_other_ops_individual_m += num_same_other_ops_individual_m;
# 		sam_line.num_correct_m_ops = num_same_m_ops_individual_m;
# 		# print num_same_m_ops_individual_m;
# 		# exit(1);

# 	# precision_m = (float(total_num_same_m_ops_individual_m) / float(ref_num_m_ops_mapped_reads)) * 100.0;
# 	precision_m = (float(total_num_same_m_ops_individual_m) / float(query_num_m_ops_mapped_reads)) * 100.0;
# 	precision_i = (float(total_num_same_i_ops_individual_m) / float(ref_num_i_ops_mapped_reads)) * 100.0;
# 	precision_d = (float(total_num_same_d_ops_individual_m) / float(ref_num_d_ops_mapped_reads)) * 100.0;

# 	recall_m = (float(total_num_same_m_ops_individual_m) / float(ref_num_m_ops_all_reads)) * 100.0;
# 	recall_i = (float(total_num_same_i_ops_individual_m) / float(ref_num_i_ops_all_reads)) * 100.0;
# 	recall_d = (float(total_num_same_d_ops_individual_m) / float(ref_num_d_ops_all_reads)) * 100.0;

# 	# if (('LAST' in sam_basename.upper()) == True and precision_m < 95.0):
# 	# 	print 'Tu sam!!!';
# 	# 	print sam_basename;
# 	# 	print 'total_num_same_m_ops_individual_m = %d' % total_num_same_m_ops_individual_m;
# 	# 	print 'ref_num_m_ops_mapped_reads = %d' % ref_num_m_ops_mapped_reads;
# 	# 	print '';
# 	# 	print '';
# 	# 	print '';

# 	if (out_summary_prefix != ''):
# 		fp_out.write('percent_correct_m\tnum_correct_m\tnum_m_ops_in_reference\n');
# 		fp_out.write('%.2f\t%.2f\t%.2f\t%.2f\n' % (precision_m, recall_m, total_num_same_m_ops_individual_m, ref_num_m_ops_mapped_reads));
# 		fp_out.close();
	
# 	sys.stderr.write('\n');

# 	return [precision_m, recall_m, total_num_same_m_ops_individual_m, ref_num_m_ops_mapped_reads, ref_num_m_ops_all_reads];



def CalculateStats(num_unique_references, sam_lines, sam_basename, sam_execution_stats, sam_modified_time, intermediate_results_path=''):
	# Sort lines by three conditions: first by name which will help us to find duplicates,
	# then by those that have the right orientation and reference name, so we can filter
	# out those that are not good, and third by the minimum distance from the reference
#	sorted_lines_by_name = sorted(sam_lines, key=lambda sam_line: (sam_line.qname, (1 - sam_line.is_correct_ref_and_orient), sam_line.min_distance));
#	sorted_lines_by_name = sorted(sam_lines, key=lambda sam_line: (sam_line.qname, (1 - sam_line.is_correct_ref_and_orient), (-sam_line.chosen_quality)));
	sorted_lines_by_name = sorted(sam_lines, key=lambda sam_line: (sam_line.qname, (-sam_line.chosen_quality)));
	# OVO SAM BIO ZAKOMENTIRAO U ZADNJOJ VERZIJI!! utility_sam.WriteSamLines(sorted_lines_by_name, ('temp/new_output-%s.txt' % sam_basename));
	if (intermediate_results_path != ''):
		utility_sam.WriteSamLines(sorted_lines_by_name, ('%s/new_output-sorted-%s_%s.txt' % (intermediate_results_path, sam_basename, SCRIPT_VERSION)));
	
	# Filter unique SAM lines, where uniqueness is defined by the qname parameter.
	# If there is more than one alignment with the same name, pick only the first one
	# because we have already sorted them both by correctness (compared to the orientation
	# and the name of the reference sequence), and by distance from the reference.
	unique_lines = [];
	i = 0;
	previous_line = None;
	num_unambiguous_reads = 0;
	num_unmapped_alignments = 0;
	num_unmapped_reads = 0;
	num_unambiguous_mapped_reads = 0;
	is_read_mapped = 0;
	count = 0;
	while i < len(sorted_lines_by_name):
		sam_line = sorted_lines_by_name[i];
		
		if (sam_line.IsMapped() == False):
			num_unmapped_alignments += 1;
		else:
			is_read_mapped = 1;
		
		if (previous_line == None or sam_line.qname != previous_line.qname):
			if (count == 1):
				num_unambiguous_reads += 1;
				num_unambiguous_mapped_reads += is_read_mapped;
			unique_lines.append(sam_line);
			sam_line.is_duplicate = False;
			num_unmapped_reads += 1 if (sam_line.IsMapped() == False) else 0;
			count = 1;
			is_read_mapped = 0;
		elif (sam_line.qname == previous_line.qname):
			count += 1;
			sam_line.is_duplicate = True;
		previous_line = sam_line;
		i += 1;
	if (count == 1):
		num_unambiguous_reads += 1;
		num_unambiguous_mapped_reads += is_read_mapped;
	
	#OVO SAM BIO ZAKOMENTIRAO U ZADNJOJ VERZIJI!!
	if (intermediate_results_path != ''):
		utility_sam.WriteSamLines(unique_lines, ('%s/new_output-unique-%s.txt' % (intermediate_results_path, sam_basename)));
	
	[true_positive, false_positive, not_mapped] = utility_sam.GetBasicStats(unique_lines, allowed_distance=52);
	total_mapped = true_positive + false_positive;
	total_not_mapped = not_mapped + num_unique_references - total_mapped;
	summary_line = '';
	summary_line += '(recall)    TP / num_unique_references = %d (%.2f%%)\n' % (true_positive, (0.0 if (num_unique_references == 0) else ((float(true_positive) / float(num_unique_references))*100.0)) );
	summary_line += '            FP / num_unique_references = %d (%.2f%%)\n' % (false_positive, (0.0 if (num_unique_references == 0) else ((float(false_positive) / float(num_unique_references))*100.0)) );
	summary_line += '(precision) TP / total_uniquely_mapped) = %d (%.2f%%)\n' % (true_positive, (0.0 if (total_mapped == 0) else ((float(true_positive) / float(total_mapped))*100.0)) );
	summary_line += '            FP / total_uniquely_mapped = %d (%.2f%%)\n' % (false_positive, (0.0 if (total_mapped == 0) else ((float(false_positive) / float(total_mapped))*100.0)) );
	summary_line += 'total_mapped = %d (%.2f%%)\n' % (total_mapped, (0.0 if (num_unique_references == 0) else (float(total_mapped) / float(num_unique_references))*100.0) );
	summary_line += 'not_mapped = %d (%.2f%%)\n' % (total_not_mapped, (0.0 if (num_unique_references == 0) else (float(total_not_mapped) / float(num_unique_references))*100.0) );
	summary_line += 'num_alignments_in_sam = %d\n' % (len(sam_lines));
	summary_line += 'num_mapped_alignments_in_sam = %d\n' % (len(sam_lines) - num_unmapped_alignments);
	summary_line += 'num_unmapped_alignments_in_sam = %d\n' % (num_unmapped_alignments);
	summary_line += 'num_reads_in_sam = %d\n' % (len(unique_lines));
	summary_line += 'num_unambiguous_reads = %d\n' % (num_unambiguous_reads);
	
	distance_limits = range(0, 60);
	#[distance_histogram_x, distance_histogram_y] = utility_sam.GetDistanceHistogramStats(unique_lines, distance_limits);
	[distance_histogram_x, distance_histogram_y] = utility_sam.GetDistanceHistogramStatsScaleDuplicates(unique_lines, distance_limits, scale_by_num_occurances=False);
	distance_histogram = [distance_histogram_x, distance_histogram_y];
	
	# Correctness is the percentage of correctly mapped alignments to within the given distance. The percentage is with respect to only the number of mapped alignments.
#	distance_accuracy_y = [((float(value) / float((len(sam_lines) - num_unmapped_alignments)))*100.0) for value in distance_histogram_y];	# This is from v3 of this script.
	# distance_accuracy_y = [((float(value) / float(len(unique_lines) - num_unmapped_reads))*100.0) for value in distance_histogram_y];		# This is new, from v4 of the script.
	distance_accuracy_y = [((float(value) / float(len(unique_lines) - num_unmapped_reads))*100.0) if (len(unique_lines) - num_unmapped_reads) else 0.0 for value in distance_histogram_y];		# This is new, from v4 of the script.
	distance_accuracy = [distance_histogram_x, distance_accuracy_y];

	# Correctness with unmapped is the percentage of correctly mapped alignments with respect to the total number of reported alignments (including those marked not mapped). These include all reported alignments + number of reads that were not reported in the SAM file.
#	max_number_of_alignments_with_unmapped = len(sam_lines) + (num_unique_references - len(unique_lines));	# This is from v3 of this script.
	max_number_of_alignments_with_unmapped = num_unique_references;											# This is new, from v4 of the script.
	distance_recall_y = [((float(value) / float(max_number_of_alignments_with_unmapped))*100.0) for value in distance_histogram_y];
	distance_recall = [distance_histogram_x, distance_recall_y];

	summary_line += '-\n';
	summary_line += 'distance_histogram[0]               = %d\n' % distance_histogram_y[0];
	summary_line += 'distance_histogram[10]              = %d\n' % distance_histogram_y[10];
	summary_line += 'distance_histogram[49]              = %d\n' % distance_histogram_y[49];

	summary_line += '-\n';
	summary_line += 'distance_accuracy[0]             = %.2f%%\n' % distance_accuracy_y[0];
	summary_line += 'distance_accuracy[10]            = %.2f%%\n' % distance_accuracy_y[10];
	summary_line += 'distance_accuracy[49]            = %.2f%%\n' % distance_accuracy_y[49];

	summary_line += '-\n';
	summary_line += 'distance_recall[0]  = %.2f%%\n' % distance_recall_y[0];
	summary_line += 'distance_recall[10] = %.2f%%\n' % distance_recall_y[10];
	summary_line += 'distance_recall[49] = %.2f%%\n' % distance_recall_y[49];
	
	#print distance_histogram;
	#print distance_histogram_percentage;
	
	mapq_histograms = [];
	mapq_histograms_percentage = [];
	distance_limits_for_mapq = [1, 50, 100];
	mapq_limits = range(0, 255);
	for distance_limit_for_mapq in distance_limits_for_mapq:
		[mapq_distance, mapq_histogram_x, mapq_histogram_y] = utility_sam.GetMapqHistogramStats(unique_lines, mapq_limits, distance_limit_for_mapq);
		mapq_histogram_percentage_y = [((float(mapq_num_correct) / float(num_unique_references))*100.0) for mapq_num_correct in mapq_histogram_y];
		mapq_histograms.append([mapq_distance, mapq_histogram_x, mapq_histogram_y]);
		mapq_histograms_percentage.append([mapq_distance, mapq_histogram_x, mapq_histogram_percentage_y]);

	summary_line += 'Execution stats:\n%s\n' % sam_execution_stats;
	summary_line += 'Modified time: %s\n' % sam_modified_time;
	
	# print summary_line;
	
	return [total_mapped, distance_histogram, distance_accuracy, distance_recall, mapq_histograms, mapq_histograms_percentage, summary_line];

def PlotROC(roc_curves, labels, xlabel_title, ylabel_title, figure_title, dataset_description, out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		fig = plt.figure();
		ax1 = plt.subplot(111);
		
		i = 0;
		while i < len(roc_curves):
			roc = roc_curves[i];
			[fpr, tpr] = roc;
			ax1.plot(fpr, tpr, label=labels[i][0]);
			i += 1;
		
		#ax1.set_xlabel('False positive rate (FPR) [%]');
		#ax1.set_ylabel('True positive rate (TPR) [%]');
		#ax1.set_xlabel('Fraction of reads [%]');
		#ax1.set_ylabel('Correctly mapped reads [%]');
		ax1.set_xlabel(xlabel_title);
		ax1.set_ylabel(ylabel_title);
		
		fontP = FontProperties()
		fontP.set_size('small')
		ax1.grid();
		
		if (len(labels) > 0):
			ax1.legend(prop=fontP, loc='lower right');
		
		#title_string = 'ROC curve of mapping position accuracy\n%s' % dataset_description;
		title_string = '%s\n%s' % (figure_title, dataset_description);
		ax1.set_title(title_string);
		
		font = {'family' : 'sans-serif',
			'weight' : 'normal',
			'size'   : 10}

		plt.rc('font', **font);
		
		#x1,x2,y1,y2 = plt.axis();
		plt.axis((-0.05, 1.05, -0.05, 1.05));
		
		if (out_png_path != ''):
			plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);

def PlotHistogramsDistance(histograms, labels, dataset_description, out_png_path):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		fig = plt.figure();
		ax1 = plt.subplot(111);
		
		#for histogram in histograms:
		i = 0;
		while i < len(histograms):
			histogram = histograms[i];
			[xvalues, yvalues] = histogram;
			xvalues = xvalues;
			yvalues = yvalues;
			ax1.plot(xvalues, yvalues, label=labels[i][0]);
			i += 1;
		
		ax1.set_xlabel('Allowed distance [bp]');
		ax1.set_ylabel('Correctly mapped reads [%]');
		
		fontP = FontProperties()
		fontP.set_size('small')
		ax1.grid();
		if (len(labels) > 0):
			ax1.legend(prop=fontP, loc='lower right');
		
		title_string = 'Mapping accuracy to within allowed distance from expected location\n%s' % dataset_description;
		ax1.set_title(title_string);
		
		font = {'family' : 'sans-serif',
			'weight' : 'normal',
			'size'   : 10}

		plt.rc('font', **font)
		
		plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);

def PlotHistogramsMapq(histograms, labels):
	fig = None;

	if USE_MATPLOTLIB == True:
		fig = plt.figure();
		ax1 = plt.subplot(111);
		
		i = 0;
		while i < len(histograms):
			histograms2 = histograms[i];
			j = 0;
			base_line = None;
			linestyles = ['-', '--', '-.', ':', 'steps'];
			while j < len(histograms2):
				histogram = histograms2[j];
				[distance, xvalues, yvalues] = histogram;
				xvalues = range(0, len(yvalues));
				if j == 0:
					base_line, = ax1.plot(xvalues, yvalues, label=(labels[i][0] + '-' + str(distance)));
				else:
					current_line = ax1.plot(xvalues, yvalues, label=(labels[i][0] + '-' + str(distance)));
					plt.setp(current_line, color=base_line.get_color(), linestyle=linestyles[j]);
				j += 1;
			i += 1;
		
		fontP = FontProperties()
		fontP.set_size('small')
		ax1.grid();
		if (len(labels) > 0):
			ax1.legend(prop=fontP, loc='lower right');
		
		font = {'family' : 'sans-serif',
			'weight' : 'normal',
			'size'   : 10}

		plt.rc('font', **font)

def WriteCorrectBaseCount(sam_path, all_correctly_mapped_bases_titles, all_correctly_mapped_bases, all_labels, machine_name, genome_name, out_path):
	fp = None;
	
	try:
		fp = open(out_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_path;
		return;
	
	# fp.write('%s\n' % sam_path);
	# fp.write('%s\n' % out_path);

	line = '';
	line += 'x\t';
	line += '-\t';
	line += 'title\t';
	line += '\t'.join(all_correctly_mapped_bases_titles);
	fp.write(line + '\n');

	i = 0;
	while i < len(all_correctly_mapped_bases):

		mapper_label = all_labels[i][0];
		
		line = '';
		line += 'y\t';
		line += '%s\t' % mapper_label;
		line += 'values\t';
		line += '\t'.join([str(value) for value in all_correctly_mapped_bases[i]]);
		fp.write(line + '\n');
		
		i += 1;
	
	fp.close();

def WriteROC(sam_path, all_roc_curves, all_labels, machine_name, genome_name, out_path):
	fp = None;
	
	try:
		fp = open(out_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_path;
		return;
	
	fp.write('%s\n' % sam_path);
	fp.write('%s\n' % out_path);
	
	i = 0;
	while i < len(all_roc_curves):
		roc = all_roc_curves[i];
		mapper_label = all_labels[i][0];
		
		[xvalues, yvalues] = roc;
		
		line = '';
		line += 'x\t';
		line += '%s\t' % mapper_label;
		line += 'Recall\t';
		line += '\t'.join([str(xvalue) for xvalue in xvalues]);
		fp.write(line + '\n');
		
		line = '';
		line += 'y\t';
		line += '%s\t' % mapper_label;
		line += 'Precision\t';
		line += '\t'.join([str(yvalue) for yvalue in yvalues]);
		fp.write(line + '\n');
		
		i += 1;
	
	fp.close();
	
def WriteHistogramsDistance(sam_path, all_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, out_path):
	fp = None;
	
	try:
		fp = open(out_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_path;
		return;
	
	fp.write('%s\n' % sam_path);
	fp.write('%s\n' % out_path);

	i = 0;
	while i < len(all_histograms):
		histogram = all_histograms[i];
		mapper_label = all_labels[i][0];
		mapper_total_mapped = all_total_mapped[i];
		
		[xvalues, yvalues] = histogram;
		if (i == 0):
			line = '';
			line += 'x\t';
			line += '-\t';
			line += 'distance\t';
			line += '\t'.join([str(xvalue) for xvalue in xvalues]);
			fp.write(line + '\n');
		
		line = '';
		line += 'y\t';
		line += '%s\t' % mapper_label;
		line += 'read_count\t';
		line += '\t'.join([str(yvalue) for yvalue in yvalues]);
		fp.write(line + '\n');
		
		i += 1;
	
	fp.close();
	

	
def WriteHistogramsMapq(sam_path, all_histograms, all_total_mapped, all_labels, num_unique_references, machine_name, genome_name, out_path):
	fp = None;
	
	try:
		fp = open(out_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_path;
		return;
	
	fp.write('%s\n' % sam_path);
	fp.write('%s\n' % out_path);

	i = 0;
	while i < len(all_histograms):
		histograms2 = all_histograms[i];
		mapper_label = all_labels[i][0];
		mapper_total_mapped = all_total_mapped[i];
		
		j = 0;
		while j < len(histograms2):
			histogram = histograms2[j];
			[distance, xvalues, yvalues] = histogram;
			
			if (i == 0 and j == 0):
				line = '';
				line += 'x\t';
				line += '-\t';
				line += 'mapq\t';
				line += '\t'.join([str(xvalue) for xvalue in xvalues]);
				fp.write(line + '\n');
		
			line = '';
			line += 'y\t';
			line += '%s\t' % mapper_label;
			line += 'read_count_at_distance_%d\t' % distance;
			line += '\t'.join([str(yvalue) for yvalue in yvalues]);
			fp.write(line + '\n');
			
			j += 1;
		i += 1;
	
	fp.close();


def WriteAccuracies(hashed_sam_lines, out_path_prefix, write_duplicates=True):
	out_path_accuracies = out_path_prefix + '_%s.acc' % (SCRIPT_VERSION);
	
	try:
		fp = open(out_path_accuracies, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_path_accuracies;
		return;
	
	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: (sam_line.IsMapped(), sam_line.min_distance));
	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: (sam_line.IsMapped(), len(sam_line.cigar)));
	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: (sam_line.IsMapped(), len(sam_line.cigar)));

	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.actual_ref_pos);
	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.IsReverse());
	# sorted_sam_lines = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.chosen_quality);
	unique_sam_lines = [];
	for qname in hashed_sam_lines.keys():
		sam_line = hashed_sam_lines[qname][0];
		unique_sam_lines.append(sam_line);
	sorted_sam_lines = sorted(unique_sam_lines, reverse=True, key=lambda sam_line: (sam_line.IsMapped(), sam_line.min_distance));
	# sorted_sam_lines = sorted(unique_sam_lines, reverse=True, key=lambda sam_line: (sam_line.IsMapped(), sam_line.num_correct_m_ops));
	for sam_line in sorted_sam_lines:
		line = sam_line.FormatAccuracy();
		fp.write(line + '\n');

	# for sam_line in sorted_sam_lines:
	# 	line = sam_line.FormatAccuracy();
	# 	if (write_duplicates == False and sam_line.is_duplicate == True):
	# 		continue;
	# 	fp.write(line + '\n');
	
	fp.close();

def VerboseSumScores(sam_filename, scores, rmsd, modified_time, out_path_prefix, execution_time, fp_sum):
	roc_fp = scores[0];
	roc_tp = scores[1];
	correct = scores[2];
	wrong = scores[3];
	not_counted = scores[4];
	not_mapped = scores[5];
	num_all_alignments = scores[6];
	num_correct_all = scores[7];
	num_wrong_all = scores[8];
	num_references = scores[9];
	
	total_mapped = correct + wrong;
	
	fp_sum.write('Statistics for file "%s":\n' % sam_filename);
	fp_sum.write('Summary file path: %s\n' % ((out_path_prefix + '.sum')));
	fp_sum.write('Accuracy file path: %s\n' % ((out_path_prefix + '.acc')));
	fp_sum.write('Modified time: ' + str(modified_time));
	fp_sum.write('\n');
	fp_sum.write('Correct: %d (%.2f%%)\n' % (correct, (100.0*float(correct)/float(correct + wrong)) if (correct + wrong) > 0 else float('NaN')));
	fp_sum.write('Wrong: %d (%.2f%%)\n' % (wrong, (100.0*float(wrong)/float(correct + wrong)) if (correct + wrong) > 0 else float('NaN')));
	fp_sum.write('Total mapped: %d (%.2f%% of aligner\'s SAM entries, %.2f%% of all reference alignments)\n' % (total_mapped, (100.0*float(total_mapped)/float(total_mapped + not_mapped)) if (total_mapped + not_mapped) > 0 else float('NaN'), (100.0*float(total_mapped)/float(num_references)) if (num_references) > 0 else float('NaN')));
	fp_sum.write('Not mapped: %d (%.2f%% of SAM entries)\n' % (not_mapped, (100.0*float(not_mapped)/float(total_mapped + not_mapped)) if (total_mapped + not_mapped) > 0 else float('NaN')));
	fp_sum.write('Not counted: %d\n' % not_counted);
	fp_sum.write('Number of all alignments in the tested SAM file: %s\n' % num_all_alignments);
	fp_sum.write('Number of all correct alignments (even duplicates): %d (%.2f%%)\n' % (num_correct_all, (100.0*float(num_correct_all)/float(num_correct_all + num_wrong_all)) if (num_correct_all + num_wrong_all) > 0 else float('NaN')));
	fp_sum.write('Number of all wrong alignments (even duplicates): %d (%.2f%%)\n' % (num_wrong_all, (100.0*float(num_wrong_all)/float(num_correct_all + num_wrong_all)) if (num_correct_all + num_wrong_all) > 0 else float('NaN')));
	fp_sum.write('Number of sequences in reference SAM: %d\n' % num_references);
	fp_sum.write('Percentage of mapped sequences compared to number of seqs in ref. SAM: %.2f%%\n' % ((100.0*float(total_mapped)/float(num_references)) if (num_references) > 0 else float('NaN')));
	fp_sum.write('Percentage of all alignments compared to number of seqs in ref. SAM: %.2f%%\n' % ((100.0*float(num_correct_all + num_wrong_all)/float(num_references)) if (num_correct_all + num_wrong_all) > 0 else float('NaN')));
	fp_sum.write('RMSD: %f\n' % rmsd);
	fp_sum.write('Execution stats:\n%s\n' % execution_time);

def VerboseRoc(scores, fp_roc):
	roc_fp = scores[0];
	roc_tp = scores[1];
	num_references = scores[9];
	
	roc_x = roc_fp + [num_references];
	roc_y = roc_tp + [num_references];

	i = 0;
	fp_roc.write('num_references: %d\n' % num_references);
	while i < len(roc_x):
		fp_roc.write('%d\t%d\n' % (roc_x[i], roc_y[i]));
		i += 1;

def WriteSummary(sam_filenames, summary_lines, out_summary_path):
	try:
		fp_sum = open(out_summary_path, 'w');
	except IOError:
		print 'ERROR: Could not open file(s) for writing!';
		return;
	
	i = 0;
	while i < len(sam_filenames):
		sam_filename = sam_filenames[i];
		summary_line = summary_lines[i];
		sam_timestamp = str(os.path.getmtime(sam_filename));
		
		fp_sum.write('SAM file: %s\n' % sam_filename);
		fp_sum.write('SAM timestamp: %s\n' % sam_timestamp);
		fp_sum.write(summary_line);
		fp_sum.write('\n');
		i += 1;

	fp_sum.close();



def verbose_usage_and_exit():
	sys.stderr.write('Compares two sam files and compares alignment positions and accuracy of mapped bases.\n');
	sys.stderr.write('\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s <alignments.sam> <reference.sam> count_bases [<output_path>]' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\tcount_bases - Can be "true" or "false". If true, base alignment accuracy will be performed as well (time consuming).\n');
	sys.stderr.write('\n');

	exit(0);

def main():
	if (len(sys.argv) < 4 or len(sys.argv) > 5):
		verbose_usage_and_exit();

	sam_file = sys.argv[1];
	reference_sam = sys.argv[2];
	count_correct_bases = True if (sys.argv[3].lower() == 'true') else False;
	output_path = os.path.dirname(sam_file);
	if (len(sys.argv) == 5):
		output_path = sys.argv[4];

	timestamp = time.strftime("%Y_%m_%d/%H%M%S");

	if (len(output_path) == 0):
		output_path = './';

	EvaluateAlignments(reference_sam, [sam_file], timestamp, output_path, count_correct_bases=count_correct_bases);



if __name__ == "__main__":
	main();
