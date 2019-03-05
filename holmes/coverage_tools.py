
from subprocess import call, Popen, PIPE

import argparse

import numpy

import sys

import shutil

import csv

import copy

import os

import itertools

from itertools import islice

import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import pylab as pl
import math
import pysam
from os import listdir
from os.path import isfile, join

def average(lst):
	return sum(lst) / float(len(lst))

def get_region_coverage_info(bam_file_path, start_pos, stop_pos, chr_name, threshold=16):

	if chr_name[:3] == 'chr':
		chr_name = chr_name[3:]

	samfile = pysam.AlignmentFile(bam_file_path, "rb")

	iter = samfile.pileup(chr_name, int(start_pos), int(stop_pos))

	data_mapped = {}
	x_data = []
	y_data = []
	seq_length = stop_pos - start_pos + 1
	coverage_total = 0
	bases_over_coverage_threshold = 0

	for x in iter:
		data_mapped[x.reference_pos] = x.nsegments

	for i in range(start_pos, stop_pos):
		if i in data_mapped.keys():
			x_data.append(i)
			y_data.append(data_mapped[i])
			coverage_total += float(data_mapped[i])

			if float(data_mapped[i]) >= threshold:
				bases_over_coverage_threshold += 1



		else:
			x_data.append(i)
			y_data.append(0)



	a_matrix = {}
	a_matrix['position'] = x_data
	a_matrix['coverage'] = y_data
	a_matrix['average_coverage'] = round((coverage_total / seq_length), 0)
	a_matrix['coverage_over_threshold'] = round((100 / float(seq_length)) * float(bases_over_coverage_threshold), 0)

	return a_matrix


def parse_annotation_file(in_file_path):

	file_obj = open(in_file_path, 'r')

	anno_obj = {}

	if in_file_path.split('.')[-1] == 'bed':
		print('Parsing bed file')
		for line in file_obj:
			line_list = line.split('\t')

			gene_name = line_list[3].split('.')[0]
			# Catch for only one exon
			if len(line_list[3].split('.')) == 1:
				exon_number = '0'
			else:
				exon_number = line_list[3].split('.')[1].split('_')[2]

			gene_chrom = line_list[0]
			gene_start = line_list[1]
			gene_stop = line_list[2]

			if gene_name in anno_obj.keys():
				anno_obj[gene_name]['exons'][exon_number] = {'exon_start': gene_start}
				anno_obj[gene_name]['exons'][exon_number]['exon_stop'] = gene_stop
				anno_obj[gene_name]['gene_chrom'] = gene_chrom

			else:
				anno_obj[gene_name] = {'exons': {exon_number: {'exon_start': gene_start, 'exon_stop': gene_stop}}, 'gene_chrom': gene_chrom}

	return anno_obj


def calc_all_gene_coverage_info(anno_file_path, bam_file_path, threshold=16):
	'''
	This may be bugged because of different chromosomes
	:param anno_file_path:
	:param bam_file_path:
	:param threshold:
	:return:
	'''
	cov_info = {}

	gene_anno = parse_annotation_file(anno_file_path)

	for a_gene in gene_anno.keys():

		cov_info[a_gene] = {}

		for an_exon in gene_anno[a_gene]['exons'].keys():

			exon_coverage = get_region_coverage_info(bam_file_path,
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_start']),
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_stop']),
													 gene_anno[a_gene]['gene_chrom'], threshold=threshold)
			cov_info[a_gene][an_exon] = exon_coverage
			cov_info[a_gene][an_exon]['exon_start'] = int(gene_anno[a_gene]['exons'][an_exon]['exon_start']),
			cov_info[a_gene][an_exon]['exon_stop'] = int(gene_anno[a_gene]['exons'][an_exon]['exon_stop']),

	return cov_info


def calc_gene_coverage_info(anno_file_path, a_gene, bam_file_dir, threshold=16):
	'''
	This may be bugged because of different chromosomes

	'''
	cov_info = {}

	bam_file_list = [f for f in listdir(bam_file_dir) if isfile(join(bam_file_dir, f)) and f[-3:] == 'bam']

	print(bam_file_list)

	gene_anno = parse_annotation_file(anno_file_path)

	for a_bam_file in bam_file_list:

		sample_name = a_bam_file.split('.')[0]
		print('Processing ' + sample_name)

		cov_info[sample_name] = {}
		bam_file_path = bam_file_dir + a_bam_file

		for an_exon in gene_anno[a_gene]['exons'].keys():

			exon_coverage = get_region_coverage_info(bam_file_path,
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_start']),
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_stop']),
													 gene_anno[a_gene]['gene_chrom'], threshold=threshold)
			cov_info[sample_name][an_exon] = exon_coverage
			cov_info[sample_name][an_exon]['exon_start'] = int(gene_anno[a_gene]['exons'][an_exon]['exon_start']),
			cov_info[sample_name][an_exon]['exon_stop'] = int(gene_anno[a_gene]['exons'][an_exon]['exon_stop']),

	return cov_info


def calc_average_gene_coverage_info(anno_file_path, a_gene, bam_file_dir, threshold=16):
	'''
	This may be bugged because of different chromosomes

	'''
	cov_info = {}

	bam_file_list = [f for f in listdir(bam_file_dir) if isfile(join(bam_file_dir, f)) and f[-3:] == 'bam']

	print(bam_file_list)

	gene_anno = parse_annotation_file(anno_file_path)

	for a_bam_file in bam_file_list:

		sample_name = a_bam_file.split('.')[0]
		print('Processing ' + sample_name)

		cov_info[sample_name] = {}
		bam_file_path = bam_file_dir + a_bam_file

		total_base_coverage_list = []

		for an_exon in gene_anno[a_gene]['exons'].keys():

			exon_coverage = get_region_coverage_info(bam_file_path,
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_start']),
													 int(gene_anno[a_gene]['exons'][an_exon]['exon_stop']),
													 gene_anno[a_gene]['gene_chrom'], threshold=threshold)

			total_base_coverage_list += exon_coverage['coverage']

		print(total_base_coverage_list)

		bases_over_threshold = 0
		for a_base in total_base_coverage_list:
			if int(a_base) > threshold:
				bases_over_threshold += 1

		print(len(total_base_coverage_list))
		print(bases_over_threshold)

		percentage_of_gene_covered_over_threshold = 100.0 / len(total_base_coverage_list) * bases_over_threshold

		print(percentage_of_gene_covered_over_threshold)
		print('average_cov')
		print(average(total_base_coverage_list))

		cov_info[sample_name]['average_coverage'] = average(total_base_coverage_list)
		cov_info[sample_name]['percentage_over_cov_thres'] = percentage_of_gene_covered_over_threshold
		cov_info[sample_name]['percentage_over_cov_thres_perc'] = cov_info[sample_name]['percentage_over_cov_thres'] / 100.0


	return cov_info


def plot_coverage_figure_for_gene(cov_result, a_gene):
	'''
	Ok, we need some standards here.
	We need the annotation positions
	We need to a coverage matrix, with positions related to depth, per sample.
	:param cov_result:
	:param a_gene:
	:return:
	'''
	exon_yposition = -10
	coverage_array = []
	positions_array = []
	exon_position_list = []

	for exon in sorted(cov_result[a_gene].keys()):
		coverage_array = coverage_array + cov_result[a_gene][exon]['coverage']
		positions_array = positions_array + cov_result[a_gene][exon]['position']
		exon_position_list.append([[cov_result[a_gene][exon]['exon_start'], cov_result[a_gene][exon]['exon_stop']],[exon_yposition, exon_yposition]])

	list1, list2 = zip(*sorted(zip(positions_array, coverage_array)))

	print(list1[-1] - list1[0])
	print(len(list1))

	full_coverage_array = []
	full_positions_array = []
	'''
	for an_x in range(cov_result[a_gene][],11):

	if an_x in positions_array:
		full_positions_array
	
	'''
	df = pd.DataFrame({'x': list1, 'y1': list2})

	plt.figure(dpi=200, figsize=[16.4, 14.8])

	plt.plot('x', 'y1', data=df, color='skyblue')


	# Draw exons and what not
	for an_exon in exon_position_list:
		print(an_exon)
		print(an_exon[0])
		plt.plot(an_exon[0], an_exon[1], 'k-', lw=2, color='red')



	plt.savefig('tesetPlot')
	plt.close()

	print(df)

	return 'Done'


def main(bam_file_path):

	print("Starting analysis")

	start_pos = 32174044
	stop_pos = 32174174
	chr_name = '22'

	anno_path = '/Volumes/External/CIDRI_Data/quick_jobs/Alina_Esterhuizen_20-11-2018/test_1_renamed.bed'

	calc_gene_coverage_info(anno_path, bam_file_path)

	'''
	mapping_test = get_coverage_info(bam_file_path, start_pos, stop_pos, chr_name)

	print(mapping_test['average_coverage'])
	print(mapping_test['coverage_over_threshold'])

	fig, ax = plt.subplots()
	ax.plot(mapping_test['x_vals'], mapping_test['y_vals'])
	plt.show()
	'''




if __name__ == '__main__':

	samtoolsPath = 'samtools'
	bedtoolsPath = '/Users/panix/Dropbox/Programs/tools/bedtools2/bin/bedtools'
	annoOutName = 'testing'
	temp_folder = './'
	bam_test_file = '/Volumes/External/CIDRI_Data/quick_jobs/Alina_Esterhuizen_20-11-2018/bams2/EE134_1HAR.unique.sort.bam'

	refGenome = '/Volumes/External/Genomes/human/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa'

	main(bam_test_file)





