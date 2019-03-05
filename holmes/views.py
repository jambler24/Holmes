from django.http import HttpResponse
from django.shortcuts import render
from django.conf import settings
from django.core.files.storage import FileSystemStorage
import json, os
import networkx as nx
from networkx.readwrite import json_graph
import coverage_tools as cov_tools
from .processing import *
from .forms import UploadFileForm
from holmes.models import Experiment, CurrentSettings
from os import listdir
from os.path import isfile, join
from bs4 import BeautifulSoup
import requests


def home(request):
	return render(request, 'index.html')


def uploads(request):
	if request.method == 'POST' and request.FILES['myfile']:
		#myfile = request.FILES['myfile']
		#fs = FileSystemStorage()
		#filename = fs.save(myfile.name, myfile)
		#uploaded_file_url = fs.url(filename)
		uploaded_file_url = 'Temp disabled'
		file_type = 'gff'

		if file_type == 'gff':

			gff2network('/Users/panix/iCloud/programs/Holmes/holmes_core/uploads/sequence.gff3',
						'/Users/panix/iCloud/programs/Holmes/holmes_core/processed')

		file_type = 'expression'

		if file_type == 'expression':

			genome_graph = nx.read_graphml('processed/test_result_network.xml')

			expression2network('uploads/gene_exp.diff', genome_graph, 'processed')

			print("processing expression data: complete")

		file_type = 'variants'

		if file_type == 'variants':
			print("processing variant data")

			genome_graph = nx.read_graphml('processed/test_exp_network.xml')


			variants2network('uploads/S5527_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')
			variants2network('uploads/S507_aligned_to_h37rv_RG_sorted_filtered_ann_unique.vcf', genome_graph, 'processed')

			print("processing variant data: complete")


		# Convert gff to GML

		return render(request, 'uploads.html', {'uploaded_file_url': uploaded_file_url})

	return render(request, 'uploads.html')


def gene_search(request):

	if request.method == 'POST':
		# Extract gene subnet
		a_gene = request.POST['gene_name']
		genome_graph = nx.read_graphml('processed/test_var_network.xml')
		graph_net = extract_subview_genome(genome_graph, a_gene)

		# extract expression of all genes in subnet
		res_info = subnet_to_json_dict(graph_net)
		json_res_info = json.dumps(res_info)
		print(json_res_info)

		# Dump to JSON for vis
		json_dict = subnet2json(graph_net)

		# data = json.dumps({"nodes": json_data["nodes"], "links": json_data["links"]})
		data = json.dumps(json_dict)

		#print(json_data)
		#outfile = open('processed/tempGraph.json', 'w')
		#outfile.write(data)
		#outfile.close()

		# Get expression data
		expression_table_data = extract_and_format_gene_expression_table(graph_net, a_gene)

		# Format

		return render(request, 'gene_view.html', {'JSON_data': data, 'table': expression_table_data, 'gene_info': json_res_info})

	else:

		return render(request, 'search_genes.html')


def gene_view(request, gene_id):
	# Extract gene subnet
	genome_graph = nx.read_graphml('processed/test_var_network.xml')
	graph_net = extract_subview_genome(genome_graph, gene_id)

	# extract expression of all genes in subnet
	res_info = subnet_to_json_dict(graph_net)
	json_res_info = json.dumps(res_info)
	print(json_res_info)

	# Dump to JSON for vis
	json_dict = subnet2json(graph_net)

	# data = json.dumps({"nodes": json_data["nodes"], "links": json_data["links"]})
	data = json.dumps(json_dict)

	# print(json_data)
	# outfile = open('processed/tempGraph.json', 'w')
	# outfile.write(data)
	# outfile.close()

	# Get expression data
	expression_table_data = extract_and_format_gene_expression_table(graph_net, gene_id)

	# Format

	return render(request, 'gene_view.html', {'JSON_data': data, 'table': expression_table_data, 'gene_info': json_res_info, 'gene_name': gene_id})


def sub_graphs(request):

	genome_graph = nx.read_graphml('processed/test_var_network.xml')

	subnet_list = extract_subnets(genome_graph, condition='507midlogtreated_5527midlogtreated')

	request.session['subnet_data'] = subnet_list

	return render(request, 'subnet_list.html', {'table_data': subnet_list})


def variant_overview(request):

	gene_graph = nx.read_graphml('processed/test_var_network.xml')
	genome_graph = nx.read_graphml('processed/test3genome.xml')

	variant_info = get_variant_info(gene_graph, genome_graph)

	#subnet_list = extract_subnets(genome_graph, condition='507midlogtreated_5527midlogtreated')

	#request.session['subnet_data'] = subnet_list

	return render(request, 'variant_info.html', {'table_data': 'this'})


def sub_graph_detail(request):

	subnet_id = request.path.rsplit('/', 1)[-1]
	subnet_list = request.session.get('subnet_data', None)

	for subnet in subnet_list:
		if subnet_id == subnet['label']:
			displayed_subnet = subnet

	return render(request, 'subnet_view.html', {'subnet_data': displayed_subnet})


def coverage_summary(request):

	anno_path = '/Users/panix/Desktop/temp_data/annotations/test_1_renamed.bed'
	bam_files_dir = '/Users/panix/Desktop/temp_data/bams_test/'

	if request.method == 'POST':

		q_gene = request.POST['a_gene_selection']
		cov_threshold = int(request.POST['coverage_threshold'])

		# Getting gene info
		target_page = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' + q_gene
		r = requests.get(target_page)
		data = r.text
		soup = BeautifulSoup(data, 'html.parser')

		mydivs = soup.find("div")
		print(mydivs)

		sample_cov_result = cov_tools.calc_gene_coverage_info(anno_path, q_gene, bam_files_dir, threshold=cov_threshold)

		exon_cot_row = []
		cot_matrix = []

		# This needs changing
		for a_sample in sample_cov_result.keys():
			rep_sample_name = a_sample

		for an_exon in sample_cov_result[rep_sample_name].keys():
			exon_cot_row.append(int(an_exon))

		exon_sample_row = sorted(exon_cot_row)

		for a_sample in sample_cov_result.keys():
			exon_cot_row = [{'cot': a_sample, 'avc': a_sample, 'cot_perc': 0}]
			for exon in exon_sample_row:
				#print(exon)
				#print(sample_cov_result[a_sample][str(exon_sample_row)]['coverage_over_threshold'])
				#print(sample_cov_result[a_sample][str(exon_sample_row)]['average_coverage'])
				cot = sample_cov_result[a_sample][str(exon)]['coverage_over_threshold']
				avc = sample_cov_result[a_sample][str(exon)]['average_coverage']
				cot_perc = cot / 100
				res_dict = {'cot': cot, 'avc': avc, 'cot_perc': cot_perc}
				exon_cot_row.append(res_dict)
			cot_matrix.append(exon_cot_row)

		exon_sample_row = ['Sample'] + exon_sample_row

		print(exon_sample_row)

		view_info = {'gene': q_gene, 'cov_threshold': cov_threshold}

		return render(request, 'coverage_summary.html', {'JSON_data': q_gene, 'exon_cot_row': exon_sample_row, 'cot_matrix': cot_matrix, 'view_info': view_info})

	else:

		anno_obj = cov_tools.parse_annotation_file(anno_path)

		gene_list = []

		for a_gene in anno_obj.keys():
			gene_list.append(a_gene)

		return render(request, 'coverage_summary_select.html', {'gene_list':gene_list})


def coverage_summary_gene(request):

	anno_folder = '/Users/panix/Desktop/temp_data/annotations/'
	bam_files_dir = '/Users/panix/Desktop/temp_data/bams_test/'

	if request.method == 'POST':

		# q_gene = request.POST['a_gene_selection']
		cov_threshold = int(request.POST['coverage_threshold'])

		print(request.POST.keys())

		anno_file = request.POST['anno_selection']

		anno_path = anno_folder + anno_file

		anno_obj = cov_tools.parse_annotation_file(anno_path)

		gene_list = []
		table_matrix = []

		for a_gene in anno_obj.keys():
			gene_list.append(a_gene)

		#cov_tools.plot_coverage_figure_for_gene(cov_dict, 'NM_001244814')

		# Get a list of samples to iterate over

		bam_file_list = [f for f in listdir(bam_files_dir) if isfile(join(bam_files_dir, f)) and f[-3:] == 'bam']
		sample_list = []

		for filename in bam_file_list:
			sample_list.append(filename.split('.')[0])

		for a_gene in gene_list:

			gene_cot_row = [{'cot': a_gene, 'avc': a_gene, 'cot_perc': 0}]
			gene_cov_result = cov_tools.calc_average_gene_coverage_info(anno_path, a_gene, bam_files_dir, threshold=cov_threshold)

			gene_coverage_total_list = []

			print(a_gene)

			for a_sample in sample_list:

				print(a_sample)
				print(gene_cov_result[a_sample].keys())

				sample_cot_row = [{'cot': gene_cov_result[a_sample]['percentage_over_cov_thres'], 'avc': gene_cov_result[a_sample]['average_coverage'], 'cot_perc': gene_cov_result[a_sample]['percentage_over_cov_thres_perc']}]

				gene_cot_row = gene_cot_row + sample_cot_row

			table_matrix.append(gene_cot_row)

		gene_list = ['Gene'] + sample_list

		view_info = {'cov_threshold': cov_threshold, 'anno_file': anno_file}

		print(gene_list)
		print(len(table_matrix))
		print(table_matrix[0])

		return render(request, 'coverage_summary_genes.html', {'gene_list': gene_list, 'cot_matrix': table_matrix, 'view_info': view_info})

	else:

		anno_file_list = [f for f in listdir(anno_folder) if isfile(join(anno_folder, f)) and f[-3:] == 'bed']

		return render(request, 'coverage_summary_select_gene.html', {'anno_list': anno_file_list})


def view_region(request):

	return render(request, 'view_read_align_region.html', {'gene_list': ''})




