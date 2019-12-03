#! /usr/bin/python

machine_names = [];
simulated_datasets = [];

simulated_datasets.append('OxfordNanopore-pbsim-observed_last-2d-1k');	machine_names.append('nanopore');
#simulated_datasets.append('Illumina-1k-single_end');                                            machine_names.append('illumina');


genomes = [];
#genomes.append('hg19_v38-chr3');
genomes.append('escherichia_coli');
#genomes.append('hg19');


simulated_datasets_grid = [];
machine_names_grid = [];
genomes_grid = [];
error_rates = [0.05, 0.10, 0.15, 0.20, 0.25];
read_lengths = [1000, 2000, 3000, 4000, 5000];
genomes_grid.append('saccharomyces_cerevisiae');
#for read_length in read_lengths:
#	for error_rate in error_rates:
# 		simulated_datasets_grid.append('GridSearch-l%dk-e%.2f' % (read_length/1000, error_rate));				machine_names_grid.append('nanopore');

