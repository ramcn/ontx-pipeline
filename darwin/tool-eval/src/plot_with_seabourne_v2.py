#! /usr/bin/python

import os;
import sys;

try:
	import numpy as np;
except Exception, e:
	USE_MATPLOTLIB = False;
	print e;
	print 'Warning: NumPy not installed!';
	exit(0);

USE_MATPLOTLIB = True;
try:
	import matplotlib.pyplot as plt;
except Exception, e:
	USE_MATPLOTLIB = False;
	print e;
	print 'Warning: Matplotlib not installed!';
	exit(0);
	
try:
	from matplotlib.font_manager import FontProperties;
except Exception, e:
	USE_MATPLOTLIB = False;
	print e;
	print 'Matplotlib problem 2!';
	
try:
	import seaborn as sns;
except Exception, e:
	USE_MATPLOTLIB = False;
	print e;
	print 'Warning: Seaborn Python module not installed!';

HIGH_DPI_PLOT = False;
HIGH_DPI_PLOT = True;




def ReadlinesWrapper(file_path):
	try:
		fp = open(file_path, 'r');
		lines = fp.readlines();
		fp.close();
		return [line.strip() for line in lines];
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!' % file_path);
	
	return [];

def WritelineWrapper(line, file_path):
	try:
		fp = open(file_path, 'w');
		fp.write(line);
		fp.close();
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!' % file_path);



# LoadPlotFile parses a plot CSV file which contains lines, each with > 3 columns:
#	x or y - denoting the axis of the parameter
#	name of the mapper or '-' if not available
#	name of the parameter
#	parameter values
# Since for each parameter, there can and will be more than one mapper reported, two
# dicts are used to store the data:
#	1. dict has the key of the atribute name, it is used to retrieve all mappers
#	   which have this atribute reported
#	2. for each key in the atribute dict, another dict is stored where keys are
#	   mapper names, and values are the split components of the original line.
def LoadPlotFile(plot_file_path, suffix=''):
	lines = ReadlinesWrapper(plot_file_path);
	if (lines == []):
		return [];
	
	atributes = {};
	
	for line in lines:
		if (len(line) == 0):
			continue;
		
		split_line = line.split('\t');
		
		if (len(split_line) < 4):
			continue;
		
		if (len(split_line[0]) != 1 or (split_line[0] in 'xXyYzZ') == False):
			continue;
		
		axis_name = split_line[0];
		mapper_name = split_line[1];
		atribute_name = split_line[2];
		atribute_values = split_line[3:];

		if (suffix != ''):
			mapper_name = mapper_name.split('-' + suffix)[0];
		
		try:
			mappers_in_atribute = atributes[atribute_name];
			mappers_in_atribute[mapper_name] = split_line;
		except Exception, e:
			mappers_in_attribute = {};
			mappers_in_attribute[mapper_name] = split_line;
			atributes[atribute_name] = mappers_in_attribute;
	
	return atributes;

def GetAtributesForPlot(atributes, x_atribute_name, y_atribute_name):
	#print atributes;
	x_value_dict = atributes[x_atribute_name];
	if (len(x_value_dict) != 1):
		print 'ERROR: X values are not selected properly! More than one occurance of parameter in the input file!';
		return [];
	x_value = x_value_dict.values()[0][3:];
		
	y_value_dict = atributes[y_atribute_name];
	y_values = [];
	labels = [];
	for key in sorted(y_value_dict.keys()):
		labels.append(key);
		y_values.append(y_value_dict[key][3:]);
	
	return [x_value, y_values, labels];

def GetAtributesForROC(atributes, x_atribute_name, y_atribute_name):
	x_value_dict = atributes[x_atribute_name];
	x_values = [];
	labels = [];
	for key in sorted(x_value_dict.keys()):
		labels.append(key);
		x_values.append(x_value_dict[key][3:]);
	
	y_value_dict = atributes[y_atribute_name];
	y_values = [];
	for key in labels:
		y_values.append(y_value_dict[key][3:]);
	
	return [x_values, y_values, labels];



#def AutoLabel(ax, rects):
    #for rect in rects:
        #h = rect.get_height()
        ##ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                ##ha='center', va='bottom')
        #ax.text(rect.get_x()+rect.get_width()/2., 1.01*h, '%d'%int(h),
                #ha='center', va='bottom')

def GetMapperName(label):
	mapper_name = label;
	
	# mapper_name = label.split('-')[0];
	# try:
	# 	mapper_name = mapper_name_lookup[mapper_name];
	# except:
	# 	mapper_name = mapper_name[0].upper() + mapper_name[1:];
	
	return mapper_name;
	
	#try:
		#mapper_name = mapper_name_lookup[label];
	#except:
		#split_label = label.split('-');
		#mapper_name = ' '.join(split_label[0:-1]);
		#mapper_name = mapper_name[0].upper() + mapper_name[1:].lower();
		#if (split_label[-1].startswith('params_')):
			#mapper_name += ' (%s)' % ((split_label[-1].split('params_')[-1]));
	
	#return mapper_name;
	
def PlotLines(x_value, y_values, labels, x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		plt.clf();
		#sns.set_style("whitegrid");
		sns.set_style("darkgrid");
		sns.set_style("white")
		sns.set_style("ticks");

		x_min_index = 0;
		x_max_index = len(x_value);
		
		if (x_min != ''):
			i = 0;
			while i < len(x_value):
				if (float(x_value[i]) >= float(x_min)):
					x_min_index = i;
					break;
				i += 1;
			
		if (x_max != ''):
			i = len(x_value) - 1;
			while i >= 0:
				if (float(x_value[i]) < float(x_max)):
					x_max_index = i + 1;
					break;
				i -= 1;

		i = 0;
		while i < len(y_values):
			mapper_name = GetMapperName(labels[i]);

			# linestyle or ls 	[ '-' | '--' | '-.' | ':' | 'steps' | ...]
			if (i < 5):
				plt.plot(x_value[x_min_index:x_max_index], y_values[i][x_min_index:x_max_index], label=mapper_name);
			elif (i < 10):
				plt.plot(x_value[x_min_index:x_max_index], y_values[i][x_min_index:x_max_index], '--', label=mapper_name);
			else:
				plt.plot(x_value[x_min_index:x_max_index], y_values[i][x_min_index:x_max_index], '-.', label=mapper_name);
			i += 1;
			
		plt.grid();
		lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
		
		plt.xlabel(x_title);
		plt.ylabel(y_title);
		#plt.title(title);

		# plt.text(0.5, 1.08, title,
		# 	horizontalalignment='center',
		# 	fontsize=12,
		# 	transform = plt.gca().transAxes)
		
		#sns.despine();
		sns.despine(offset=10, trim=True);
		
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

def AutoLabel(ax, rects, prefix='', suffix=''):
    for rect in rects:
        h = rect.get_height()
        #ax.text(rect.get_x()+rect.get_width()/2., 1.01*h, '%s%d%s' % (prefix, int(h), suffix),
                #ha='center', va='bottom', fontsize=5)
        #ax.text(rect.get_x()+rect.get_width()/2., 1.00*h + (-25), '%s' % (prefix),
                #ha='center', va='bottom', fontsize=5)
        #ax.text(rect.get_x()+rect.get_width()/2., 1.00*h + 100, '%d%s' % (int(h), suffix),
                #ha='center', va='bottom', fontsize=10)
	
        ax.text(rect.get_x()+rect.get_width()/2., 1.00*h + (0), '%s' % (prefix),
                ha='center', va='bottom', fontsize=5)
        ax.text(rect.get_x()+rect.get_width()/2., 1.00*h + 1, '%d%s' % (int(h), suffix),
                ha='center', va='bottom', fontsize=10)

def PlotBars(x_value, y_values, labels, x_param_name, x_title='X', y_title='Y', title='', out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		plt.clf();
		#sns.set_style("whitegrid");
		sns.set_style("darkgrid");
		sns.set_style("white")
		sns.set_style("ticks");
		
		multihist_width = 0.75;
		bar_width = multihist_width / len(y_values);
		center = (multihist_width - bar_width) / 2;
		#color_list = ['b', 'g', 'r', 'm', 'k', 'y'];
		
		x_index = 0;
		
		if (x_param_name == ''):
			return;

		i = 0;
		while i < len(x_value):
			if (x_value[i] == x_param_name):
				x_index = i;
				break;
			i += 1;

		x_value_num = np.arange(len(x_value));
		color_list = plt.rcParams['axes.color_cycle']
		
		i = 0;
		print 'y_values: ', y_values;
		while i < len(y_values):
			print ('y_values[%d]: ' % i), y_values[i];
			y_value_num = np.array([float(value) for value in y_values[i]]);
			print 'y_value_num: ', y_value_num;
			
			mapper_name = GetMapperName(labels[i]);
			rect = plt.bar(i, y_value_num[x_index], width=bar_width, color=color_list[i],align='center', label=mapper_name);
			AutoLabel(plt, rect, ('%' if '%' in y_title else ''));
			i += 1;
		
		plt.grid();
		#lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
		
		plt.xlabel(x_title);
		plt.ylabel(y_title);
		#plt.title(title);

		plt.text(0.5, 1.08, title,
			horizontalalignment='center',
			fontsize=12,
			transform = plt.gca().transAxes)
		
		plt.gca().set_xticks(range(len(y_values)))
		plt.gca().set_xticklabels([GetMapperName(label) for label in labels])
		
		sns.despine(offset=10, trim=True);
		
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

def PlotBarsSimple(y_values, labels, x_title='X', y_title='Y', title='', out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		plt.clf();
		#sns.set_style("whitegrid");
		sns.set_style("darkgrid");
		sns.set_style("white")
		sns.set_style("ticks");
		
		multihist_width = 0.75;
		bar_width = multihist_width / len(y_values);
		center = (multihist_width - bar_width) / 2;
		#color_list = ['b', 'g', 'r', 'm', 'k', 'y'];

		color_list = plt.rcParams['axes.color_cycle']
		
		suffix = '';
		if ('[x]' in y_title):
			suffix = 'x';
		elif ('[%]' in y_title):
			suffix = '%';
		
		i = 0;
		print 'y_values: ', y_values;
		while i < len(y_values):
			mapper_name = GetMapperName(labels[i]);
			rects = plt.bar(i, y_values[i], width=bar_width, color=color_list[i],align='center', label=mapper_name);
			#AutoLabel(plt, rect, ('%' if '%' in y_title else ''));
			
			for rect in rects:
				h = rect.get_height()
				plt.text(rect.get_x()+rect.get_width()/2., 1.01*h + 0.5, '%.2f%s' % (float(h), suffix),
					ha='center', va='bottom', fontsize=10)
			
			i += 1;
		
		plt.grid();
		lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
		
		plt.xlabel(x_title);
		plt.ylabel(y_title);
		#plt.title(title);

		plt.text(0.5, 1.08, title,
			horizontalalignment='center',
			fontsize=12,
			transform = plt.gca().transAxes)
		
		plt.gca().set_xticks(range(len(y_values)))
		plt.gca().set_xticklabels([GetMapperName(label) for label in labels])
		
		sns.despine(offset=10, trim=True);
		
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

def PlotMultiBars(x_value, y_values, labels, x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		plt.clf();
		#sns.set_style("whitegrid");
		sns.set_style("darkgrid");
		sns.set_style("white")
		sns.set_style("ticks");
		
		#multihist_width = 0.4;
		multihist_width = 0.85;
		bar_width = multihist_width / len(y_values);
		center = (multihist_width - bar_width) / 2;
		#color_list = ['b', 'g', 'r', 'm', 'k', 'y'];
		
		x_min_index = 0;
		x_max_index = len(x_value);
		
		if (x_min != ''):
			#try:
				#x_min_float = float(x_min);
				#x_max_float = float(x_max);
			i = 0;
			while i < len(x_value):
				if (float(x_value[i]) >= float(x_min)):
					x_min_index = i;
					break;
				i += 1;
			
		if (x_max != ''):
			i = len(x_value) - 1;
			while i >= 0:
				if (float(x_value[i]) < float(x_max)):
					x_max_index = i + 1;
					break;
				i -= 1;

		print 'x_min_index = ', x_min_index;
		print 'x_max_index = ', x_max_index;
		
		x_value_num = np.arange(len(x_value[x_min_index:x_max_index]));
		color_list = plt.rcParams['axes.color_cycle']
		
		min_y = 0.0;
		max_y = 0.0;
		
		#font = {'family' : 'sans-serif',
			#'weight' : 'normal',
			#'size'   : 10}
		#plt.rc('font', **font)
		
		i = 0;
		while i < len(y_values):
			y_value_num = np.array([int(value) for value in y_values[i]]);
			print 'y_value_num[x_min_index:x_max_index]: ', y_value_num[x_min_index:x_max_index];
			current_min_y = min(y_value_num[x_min_index:x_max_index]);
			current_max_y = max(y_value_num[x_min_index:x_max_index]);
			if (i == 0):
				min_y = current_min_y;
				max_y = current_max_y;
			if (current_min_y < min_y):
				min_y = current_min_y;
			if (current_max_y > max_y):
				max_y = current_max_y;
			
			mapper_name = GetMapperName(labels[i]);
			#rects = plt.bar(x_value_num[x_min_index:x_max_index] + bar_width * i - center, y_value_num[x_min_index:x_max_index], width=bar_width, color=color_list[i],align='center', label=mapper_name);
			rects = plt.bar(x_value_num + bar_width * i - center, y_value_num[x_min_index:x_max_index], width=bar_width, color=color_list[i],align='center', label=mapper_name);
			#rects = plt.bar(x_value_num + bar_width * i - center, y_value_num, width=bar_width, align='center', label=labels[i]);
			#AutoLabel(plt, rects);
			AutoLabel(plt, rects, prefix=(mapper_name + '\n'), suffix=('%' if '%' in y_title else ''));
			i += 1;
		
		plt.grid();
		lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
		
		plt.xlabel(x_title);
		plt.ylabel(y_title);
		#plt.title(title);

		plt.text(0.5, 1.08, title,
			horizontalalignment='center',
			fontsize=12,
			transform = plt.gca().transAxes)

		#plt.axis((float(0), float(len(x_value_num)), float(min_y), float(max_y)*1.15));
		
		#plt.gca().set_xticks(x_value_num[x_min_index:x_max_index] + center)
		#plt.gca().set_xticklabels(x_value[x_min_index:x_max_index])
		
		#plt.gca().set_xticks(x_value_num + center)
		#plt.gca().set_xticks(x_value_num[x_min_index:x_max_index])
		#plt.gca().set_xticklabels(x_value[x_min_index:x_max_index])
		plt.gca().set_xticks(x_value_num)
		plt.gca().set_xticklabels(x_value[x_min_index:x_max_index])
		
		print max_y;
		#plt.gca().set_yticks(range(max_y + 1))
		
		[xmin, xmax, ymin, ymax] = plt.gca().axis();
		plt.gca().axis((-0.5, xmax, ymin, ymax));
		#print title;
		#print [xmin, xmax, ymin, ymax];
		
		#sns.despine();
		sns.despine(offset=10, trim=True);
		
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);
		
		print plt.gcf().dpi;

def PlotROC(x_values, y_values, labels, x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	fig = None;
	
	if USE_MATPLOTLIB == True:
		plt.clf();
		#sns.set_style("whitegrid");
		sns.set_style("darkgrid");
		sns.set_style("white")
		sns.set_style("ticks");

		#x_min_index = 0;
		#x_max_index = len(x_value);
		
		#if (x_min != ''):
			#i = 0;
			#while i < len(x_value):
				#if (float(x_value[i]) >= float(x_min)):
					#x_min_index = i;
					#break;
				#i += 1;
			
		#if (x_max != ''):
			#i = len(x_value) - 1;
			#while i >= 0:
				#if (float(x_value[i]) < float(x_max)):
					#x_max_index = i + 1;
					#break;
				#i -= 1;
	
	
		i = 0;
		while i < len(x_values):
			recall = x_values[i];
			precision = y_values[i];
			mapper_name = GetMapperName(labels[i]);
			plt.plot(recall, precision, label=mapper_name);
			i += 1;
		
		plt.grid();
		lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
		
		plt.xlabel(x_title);
		plt.ylabel(y_title);

		# plt.text(0.5, 1.08, title,
		# 	horizontalalignment='center',
		# 	fontsize=12,
		# 	transform = plt.gca().transAxes)
		
		#plt.axis((-0.05, 1.05, -0.05, 1.05));
		plt.axis((0.0, 1.0, 0.0, 1.0));
		
		#sns.despine();
		sns.despine(offset=10, trim=True);
		
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

	else:
		print 'Cannot use MATPLOTLIB!';
		


def ParseAndPlotLines(data_path, x_atribute, y_atribute, suffix='', x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path, suffix);
	[x_value, y_values, labels] = GetAtributesForPlot(atributes, x_atribute, y_atribute);
	PlotLines(x_value, y_values, labels, x_min, x_max, x_title, y_title, title, out_png_path);

def ParseAndPlotMultiBars(data_path, x_atribute, y_atribute, suffix='', x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path, suffix);
	[x_value, y_values, labels] = GetAtributesForPlot(atributes, x_atribute, y_atribute);
	PlotMultiBars(x_value, y_values, labels, x_min, x_max, x_title, y_title, title, out_png_path);

def ParseAndPlotBars(data_path, x_atribute, y_atribute, x_value_param, suffix='', x_title='X', y_title='Y', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path, suffix);
	[x_value, y_values, labels] = GetAtributesForPlot(atributes, x_atribute, y_atribute);
	PlotBars(x_value, y_values, labels, x_value_param, x_title, y_title, title, out_png_path);

def ParseAndPlotROC(data_path, x_atribute, y_atribute, suffix='', x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path, suffix);
	[x_value, y_values, labels] = GetAtributesForROC(atributes, x_atribute, y_atribute);
	PlotROC(x_value, y_values, labels, x_min, x_max, x_title, y_title, title, out_png_path);


def PlotSimulatedResults(data_path_prefix, sam_suffix, genome_name, script_version='v0'):
	data_path_folder = os.path.dirname(data_path_prefix);
	# data_path_out_folder = '%s/analysis-final' % (data_path_folder);
	data_path_out_folder = data_path_folder;
	data_basename = os.path.basename(data_path_prefix); # (os.path.splitext(os.path.basename(data_path_prefix))[0]).split('-')[-1];
	
	if not os.path.exists(data_path_out_folder):
		print 'Creating output folders on path "%s".' % data_path_out_folder;
		os.makedirs(data_path_out_folder);
	
	results_correctness_path = '%s-distance_correct_mapped.csv' % (data_path_prefix);
	results_correctness_w_unmapped_path = '%s-distance_correct_w_unmapped.csv' % (data_path_prefix);
	results_precision_recall_path = '%s-precision_recall.csv' % (data_path_prefix);
	
	out_correctness_png_path = '%s/accuracy_%s-%s-%s-%s.png' % (data_path_out_folder, script_version, data_basename, genome_name, sam_suffix);
	out_correctness_w_unmapped_png_path = '%s/correctness_%s_w_unmapped-%s-%s-%s.png' % (data_path_out_folder, script_version, data_basename, genome_name, sam_suffix);
	out_precision_recall_png_path = '%s/precision_recall_%s-%s-%s-%s.png' % (data_path_out_folder, script_version, data_basename, genome_name, sam_suffix);
	
	ParseAndPlotLines(	results_correctness_path,
				'distance', 'read_count', suffix=sam_suffix,
				x_min='', x_max='',
				x_title='Distance from expected location [bases]', y_title='Correctness [%]', title='Mapping correctness wrt. to the total number mapped alignments\nSimulated reads from %s genome' % (genome_name),
				out_png_path=out_correctness_png_path);
	
	ParseAndPlotLines(	results_correctness_w_unmapped_path,
				'distance', 'read_count', suffix=sam_suffix,
				x_min='', x_max='',
				x_title='Distance from expected location [bases]', y_title='Correctness [%]', title='Mapping correctness wrt. to the total number of alignments (incl. unmapped)\nSimulated reads from %s genome' % (genome_name),
				out_png_path=out_correctness_w_unmapped_png_path);

	ParseAndPlotROC(	results_precision_recall_path,
				'Recall', 'Precision', suffix=sam_suffix,
				x_min='', x_max='',
				x_title='Recall', y_title='Precision', title='Precision-Recall with strict allowed distance of 0 bases\nSimulated reads from %s genome' % (genome_name),
				out_png_path=out_precision_recall_png_path);

def ParseAndPlotConsensus(data_path_coverage, data_path_counts, x_title='X', y_title='Y', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path_coverage);
	#				x_equal='20',
#				x_title='', y_title='Variant count', title='',
				#out_png_path='%s/real-seaborn-%s-threshold_20' % (output_folder, data_basename));

	[x_cov_threshold, y_snps, labels_snps] = GetAtributesForPlot(atributes, 'coverage_threshold', 'snp_count');
	[x_cov_threshold, y_insertions, labels_insertions] = GetAtributesForPlot(atributes, 'coverage_threshold', 'insertion_count');
	[x_cov_threshold, y_deletions, labels_deletions] = GetAtributesForPlot(atributes, 'coverage_threshold', 'deletion_count');
	[x_cov_threshold, y_undercovered, labels_undercovered] = GetAtributesForPlot(atributes, 'coverage_threshold', 'num_undercovered_bases');
	[x_cov_threshold, y_average_coverage, labels_average_coverage] = GetAtributesForPlot(atributes, 'coverage_threshold', 'average_coverage');
	
	print x_cov_threshold;
	selected_column = x_cov_threshold.index('20');
	
	print 'y_snps: ', y_snps;
	print ' ';
	print 'y_insertions: ', y_insertions;
	print ' ';
	print 'y_deletions: ', y_deletions;
	print ' ';
	print 'y_undercovered: ', y_undercovered;
	print ' ';
	print 'selected_column: ', selected_column;
	print ' ';
	print 'labels_snps: ', labels_snps;
	print ' ';
	
	# x values will be SNPs, insertions, deletions and count of mapped reads
	x = ['SNP', 'Insertion', 'Deletion']; # , 'Under threshold'];
	
	# y values will be 4 lists, one list for each variant type, and each list consisting of values for different mappers. All values are for a selected coverage threshold.
	y = [[] for i in range(len(y_snps))];
	i = 0;
	print len(y_snps)
	print ' ';
	while i < len(y_snps):
		print i;
		print selected_column;
		print y_snps[i]
		print y_snps[i][selected_column]
		y[i].append(y_snps[i][selected_column]);
		y[i].append(y_insertions[i][selected_column]);
		y[i].append(y_deletions[i][selected_column]);
		#y[i].append(y_undercovered[i][selected_column]);
		i += 1;
	labels = labels_snps;
	
	print 'x: ', x;
	print 'y: ', y;
	print 'labels: ', labels;
	
	order_of_mappers_name = ['graphmap-params_default', 'bwamem-params_pacbio', 'blasr-params_pacbio', 'lastal-params_q1']; # , 'lastal-params_q2'];
	order_of_mappers_index = [labels.index(label) for label in order_of_mappers_name if (label in labels)];
	print order_of_mappers_name;
	print order_of_mappers_index;
	#x = [x[i] for i in order_of_mappers_index];
	y = [y[i] for i in order_of_mappers_index];
	labels = [labels[i] for i in order_of_mappers_index];
	y_average_coverage = [y_average_coverage[i][0] for i in order_of_mappers_index];
	
	print 'x: ', x;
	print 'y: ', y;
	print 'labels: ', labels;
	
	print 'y_average_coverage: ', y_average_coverage;
	
	

	PlotMultiBars(x, y, labels, '', '', '', y_title, title, out_png_path + '-consensus_threshold_20.png');
	PlotBarsSimple([float(value) for value in y_average_coverage], labels, '', 'Average coverage [x]', title, out_png_path + '-coverage.png');
	
#def PlotBars(x_value, y_values, labels, x_param_name, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotBarsSimple(y_values, labels, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotBars(x_value, y_values, labels, x_param_name, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotMultiBars(x_value, y_values, labels, x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):

	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-read_counts-ecoliR7.3.csv';
	atributes = LoadPlotFile(data_path_counts);
	[x_counts, y_counts, labels_counts] = GetAtributesForPlot(atributes, 'counts', 'N');
	if ('percent_mapped_reads_in_sam_compared_to_reference_sam' in x_counts):
		x_counts_index = x_counts.index('percent_mapped_reads_in_sam_compared_to_reference_sam');
		labels = labels_counts;
		print x_counts;
		print y_counts;
		print x_counts_index;
		order_of_mappers_name = ['graphmap-params_default', 'bwamem-params_pacbio', 'blasr-params_pacbio', 'lastal-params_q1'];
		order_of_mappers_index = [labels.index(label) for label in order_of_mappers_name if (label in labels)];
		
		y = [y_counts[i] for i in order_of_mappers_index];
		labels = [labels[i] for i in order_of_mappers_index];
		y = [float(values[x_counts_index]) for values in y];
		
		print y;
		
		PlotBarsSimple(y, labels, '', 'Percentage of mapped reads [%]', title, out_png_path + '-count_mapped.png');
	
	#output_folder = 'plots-2.0';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotBars(	data_path,
				#'counts', 'N',
				#x_value_param='percent_mapped_reads_in_sam_compared_to_reference_sam',
				#x_title='', y_title='Percentage of reads [%]', title='',
				#out_png_path='%s/count-seaborn-%s-num_reads.png' % (output_folder, data_basename));

	
def ParseAndPlotConsensusLambda(data_path_coverage, data_path_counts, x_title='', y_title='', title='', out_png_path=''):
	atributes = LoadPlotFile(data_path_coverage);
	#				x_equal='20',
#				x_title='', y_title='Variant count', title='',
				#out_png_path='%s/real-seaborn-%s-threshold_20' % (output_folder, data_basename));

	[x_cov_threshold, y_snps, labels_snps] = GetAtributesForPlot(atributes, 'coverage_threshold', 'snp_count');
	[x_cov_threshold, y_insertions, labels_insertions] = GetAtributesForPlot(atributes, 'coverage_threshold', 'insertion_count');
	[x_cov_threshold, y_deletions, labels_deletions] = GetAtributesForPlot(atributes, 'coverage_threshold', 'deletion_count');
	[x_cov_threshold, y_undercovered, labels_undercovered] = GetAtributesForPlot(atributes, 'coverage_threshold', 'num_undercovered_bases');
	[x_cov_threshold, y_average_coverage, labels_average_coverage] = GetAtributesForPlot(atributes, 'coverage_threshold', 'average_coverage');
	
	print x_cov_threshold;
	selected_column = x_cov_threshold.index('20');
	
	print 'y_snps: ', y_snps;
	print ' ';
	print 'y_insertions: ', y_insertions;
	print ' ';
	print 'y_deletions: ', y_deletions;
	print ' ';
	print 'y_undercovered: ', y_undercovered;
	print ' ';
	print 'selected_column: ', selected_column;
	print ' ';
	print 'labels_snps: ', labels_snps;
	print ' ';
	
	# x values will be SNPs, insertions, deletions and count of mapped reads
	x = ['SNP', 'Insertion', 'Deletion']; # , 'Under threshold'];
	
	# y values will be 4 lists, one list for each variant type, and each list consisting of values for different mappers. All values are for a selected coverage threshold.
	y = [[] for i in range(len(y_snps))];
	i = 0;
	print len(y_snps)
	print ' ';
	while i < len(y_snps):
		print i;
		print selected_column;
		print y_snps[i]
		print y_snps[i][selected_column]
		y[i].append(y_snps[i][selected_column]);
		y[i].append(y_insertions[i][selected_column]);
		y[i].append(y_deletions[i][selected_column]);
		#y[i].append(y_undercovered[i][selected_column]);
		i += 1;
	labels = labels_snps;
	
	print 'x: ', x;
	print 'y: ', y;
	print 'labels: ', labels;
	
	order_of_mappers_name = ['graphmap-params_default', 'bwamem-params_pacbio', 'blasr-params_pacbio', 'lastal-params_q1']; # , 'lastal-params_q2'];
	order_of_mappers_index = [labels.index(label) for label in order_of_mappers_name if (label in labels)];
	print order_of_mappers_name;
	print order_of_mappers_index;
	#x = [x[i] for i in order_of_mappers_index];
	y = [y[i] for i in order_of_mappers_index];
	labels = [labels[i] for i in order_of_mappers_index];
	y_average_coverage = [y_average_coverage[i][0] for i in order_of_mappers_index];
	
	print 'x: ', x;
	print 'y: ', y;
	print 'labels: ', labels;
	
	print 'y_average_coverage: ', y_average_coverage;
	
	

	PlotMultiBars(x, y, labels, '', '', '', 'Variant count', title, out_png_path + '-consensus_threshold_20.png');
	PlotBarsSimple([float(value) for value in y_average_coverage], labels, '', 'Average coverage [x]', title, out_png_path + '-coverage.png');
	
#def PlotBars(x_value, y_values, labels, x_param_name, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotBarsSimple(y_values, labels, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotBars(x_value, y_values, labels, x_param_name, x_title='X', y_title='Y', title='', out_png_path=''):
#def PlotMultiBars(x_value, y_values, labels, x_min='', x_max='', x_title='X', y_title='Y', title='', out_png_path=''):

	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-read_counts-ecoliR7.3.csv';
	atributes = LoadPlotFile(data_path_counts);
	[x_counts, y_counts, labels_counts] = GetAtributesForPlot(atributes, 'counts', 'N');
	if ('percent_mapped_reads_in_sam_compared_to_reference_sam' in x_counts):
		x_counts_index = x_counts.index('percent_mapped_reads_in_sam_compared_to_reference_sam');
		labels = labels_counts;
		print x_counts;
		print y_counts;
		print x_counts_index;
		order_of_mappers_name = ['graphmap-params_default', 'bwamem-params_pacbio', 'blasr-params_pacbio', 'lastal-params_q1'];
		order_of_mappers_index = [labels.index(label) for label in order_of_mappers_name if (label in labels)];
		
		y = [y_counts[i] for i in order_of_mappers_index];
		labels = [labels[i] for i in order_of_mappers_index];
		y = [float(values[x_counts_index]) for values in y];
		
		print y;
		
		PlotBarsSimple(y, labels, '', 'Percentage of mapped reads [%]', title, out_png_path + '-count_mapped.png');
	
	#output_folder = 'plots-2.0';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotBars(	data_path,
				#'counts', 'N',
				#x_value_param='percent_mapped_reads_in_sam_compared_to_reference_sam',
				#x_title='', y_title='Percentage of reads [%]', title='',
				#out_png_path='%s/count-seaborn-%s-num_reads.png' % (output_folder, data_basename));
	

	
if __name__ == "__main__":
	output_folder = 'plots-2.0';
	if not os.path.exists(output_folder):
		print 'Creating output folders on path "%s".' % output_folder;
		os.makedirs(output_folder);
	
	#analysis_path = '/home/ivan/work/eclipse-workspace/golden-bundle/results/reads-simulated/OxfordNanopore-pbsim-40_percent/escherichia_coli';
	#data_path = analysis_path + '/total-OxfordNanopore-pbsim-40_percent-escherichia_coli_40perc-distance_correct_mapped_40perc.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotLines(	data_path,
				#'distance', 'read_count',
				#x_min='', x_max='',
				#x_title='Distance from expected location [bases]', y_title='Correctness [%]', title='Mapping correctness wrt. to the total number of alignments marked mapped\nE. Coli simulated reads with 40% error rate',
				#out_png_path='%s/correctness-%s.png' % (output_folder, data_basename));
	
	#analysis_path = '/home/ivan/work/eclipse-workspace/golden-bundle/results/reads-simulated/OxfordNanopore-pbsim-40_percent/escherichia_coli';
	#data_path = analysis_path + '/total-OxfordNanopore-pbsim-40_percent-escherichia_coli_40perc-distance_correct_w_unmapped_40perc.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotLines(	data_path,
				#'distance', 'read_count',
				#x_min='', x_max='',
				#x_title='Distance from expected location [bases]', y_title='Correctness [%]', title='Mapping correctness wrt. to the total number of alignments (incl. unmapped)\nE. Coli simulated reads with 40% error rate',
				#out_png_path='%s/correctness-%s.png' % (output_folder, data_basename));
	
	#analysis_path = '/home/ivan/work/eclipse-workspace/golden-bundle/results/reads-simulated/OxfordNanopore-pbsim-40_percent/escherichia_coli';
	#data_path = analysis_path + '/total-OxfordNanopore-pbsim-40_percent-escherichia_coli_40perc-precision_recall_40perc.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotROC(	data_path,
				#'Recall', 'Precision',
				#x_min='', x_max='',
				#x_title='Recall', y_title='Precision', title='Precision-Recall with strict allowed distance of 0 bases\nE. Coli simulated reads with 40% error rate',
				#out_png_path='%s/precision_recall-%s.png' % (output_folder, data_basename));
	
	
	
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-consensus-lambda_reads_enumerated.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'snp_count',
				#x_min='20', x_max='20',
				#x_title='Coverage threshold (read count)', y_title='SNP count', title='SNP variants\nMikheyev & Tin Lambda dataset',
				#out_png_path='%s/consensus-seaborn-%s-snp.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-consensus-lambda_reads_enumerated.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'insertion_count',
				#x_min='20', x_max='20',
				#x_title='Coverage threshold (read count)', y_title='Insertion count', title='Insertion variants\nMikheyev & Tin Lambda dataset',
				#out_png_path='%s/consensus-seaborn-%s-insertion.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-consensus-lambda_reads_enumerated.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'deletion_count',
				#x_min='20', x_max='20',
				#x_title='Coverage threshold (read count)', y_title='Deletion count', title='Deletion variants\nMikheyev & Tin Lambda dataset',
				#out_png_path='%s/consensus-seaborn-%s-deletion.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-read_counts-lambda_reads_enumerated.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotBars(	data_path,
				#'counts', 'N',
				#x_value_param='percent_mapped_reads_in_sam',
				#x_title='', y_title='Percentage of reads [%]', title='Percentage of mapped reads\nMikheyev & Tin Lambda dataset',
				#out_png_path='%s/count-seaborn-%s-num_reads.png' % (output_folder, data_basename));

	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-consensus-ecoliR7.3.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotConsensus(
				#'/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-consensus-ecoliR7.3.csv',
				#'/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-read_counts-ecoliR7.3.csv',
##				x_equal='20',
##				x_title='', y_title='Variant count', title='',
				##out_png_path='%s/real-seaborn-%s-threshold_20' % (output_folder, data_basename));
				#out_png_path='%s/real-seaborn-ecoliR7.3' % (output_folder));

	ParseAndPlotConsensusLambda(
				'/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-consensus-lambda_reads_enumerated.csv',
				'/home/ivan/work/eclipse-workspace/data/minion-review/Lambda-Mikheyev_and_Tin/alignment/lambda_reads_enumerated/analysis-final/plot-read_counts-lambda_reads_enumerated.csv',
				out_png_path='%s/real-seaborn-lambda_reads_enumerated' % (output_folder));
	
	
	

	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-consensus-ecoliR7.3.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'snp_count',
				#x_min='20', x_max='30',
				#x_title='Coverage threshold (read count)', y_title='SNP count', title='SNP variants\nE. Coli R7.3 dataset',
				#out_png_path='%s/consensus-seaborn-%s-snp.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-consensus-ecoliR7.3.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'insertion_count',
				#x_min='20', x_max='30',
				#x_title='Coverage threshold (read count)', y_title='Insertion count', title='Insertion variants\nE. Coli R7.3 dataset',
				#out_png_path='%s/consensus-seaborn-%s-insertion.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-consensus-ecoliR7.3.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotMultiBars(	data_path,
				#'coverage_threshold', 'deletion_count',
				#x_min='20', x_max='30',
				#x_title='Coverage threshold (read count)', y_title='Deletion count', title='Deletion variants\nE. Coli R7.3 dataset',
				#out_png_path='%s/consensus-seaborn-%s-deletion.png' % (output_folder, data_basename));
	#data_path = '/home/ivan/work/eclipse-workspace/data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/analysis-final/plot-read_counts-ecoliR7.3.csv';
	#data_basename = (os.path.splitext(os.path.basename(data_path))[0]).split('-')[-1];
	#ParseAndPlotBars(	data_path,
				#'counts', 'N',
				#x_value_param='percent_mapped_reads_in_sam',
				#x_title='', y_title='Percentage of reads [%]', title='Percentage of mapped reads\nE. Coli R7.3 dataset',
				#out_png_path='%s/count-seaborn-%s-num_reads.png' % (output_folder, data_basename));

	plt.show();
	