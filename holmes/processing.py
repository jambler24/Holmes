import networkx as nx
import csv, os
import django_tables2 as tables
from holmes.models import CurrentSettings, Experiment
import re
from gengraph import *
from Bio.Seq import translate
from os import listdir
from os.path import isfile, join


def subnet2json(in_graph):
	# data = json.dumps({"nodes": json_data["nodes"], "links": json_data["links"]})

	test_dict = {"nodes": [{"name":"node1","group":1}, {"name":"node2","group":2}, {"name":"node3","group":2}, {"name":"node4","group":3}], "links": [{"source":2,"target":1,"weight":1}, {"source":0,"target":2,"weight":3}]}

	return test_dict

'''

{
  "nodes":[
		{"name":"node1","group":1},
		{"name":"node2","group":2},
		{"name":"node3","group":2},
		{"name":"node4","group":3}
	],
	"links":[
		{"source":2,"target":1,"weight":1},
		{"source":0,"target":2,"weight":3}
	]
}

'''


def input_parser(file_path):
	if file_path[-3:] == ".fa" or file_path[-6:] == ".fasta":
		input_file = open(file_path, "r")
		output_list = []
		# set variables
		sequence_details = ""
		sequence = ""

		for line in input_file:
			if line[0] == ">":
				if len(sequence_details) > 0:
					sequence_details = sequence_details.strip()
					sequence = sequence.strip()
					sequence = sequence.replace("\n", "")
					gene_ID_dict = {"gene_details": sequence_details[1:], "DNA_seq": sequence}
					output_list.append(gene_ID_dict)
					sequence_details = ""
					sequence = ""
				sequence_details = line
			else:
				sequence += line

		sequence_details = sequence_details.strip()

		sequence = sequence.strip()
		sequence = sequence.replace("\n", "")

		gene_ID_dict = {"gene_details": sequence_details[1:], "DNA_seq": sequence}

		output_list.append(gene_ID_dict)

		return output_list

	if file_path[-4:] == ".bim":
		input_file = open(file_path, "r")
		return input_file

	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result = numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".ssv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=';')
		data_matrix = list(data_table)
		result = numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".txt":
		list_of_dicts = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts

	if file_path[-4:] == ".vcf":
		list_of_dicts = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if line[0:2] != "##" and line[0] == "#":
				vcf_headder_line = line.split('\t')
				vcf_headder_line[0] = vcf_headder_line[0][1:]
				vcf_headder_line[-1] = vcf_headder_line[-1].strip()

			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {vcf_headder_line[0]: entries[0], vcf_headder_line[1]: entries[1],
							  vcf_headder_line[2]: entries[2], vcf_headder_line[3]: entries[3],
							  vcf_headder_line[4]: entries[4], vcf_headder_line[5]: entries[5],
							  vcf_headder_line[6]: entries[6], vcf_headder_line[7]: entries[7],
							  vcf_headder_line[8]: entries[8], vcf_headder_line[9]: entries[9].strip(),
							  'ORIGIN': entry_label}

				list_of_dicts.append(entry_dict)
		return list_of_dicts

	if file_path[-5:] == ".diff":
		list_of_dicts = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts

	if file_path[-5:] == "kbone":
		backb_listOlists = []
		in_file = open(file_path, 'r')
		for line in in_file:
			curr_list = []
			line = line.strip()
			curr_list = line.split('\t')
			backb_listOlists.append(curr_list)
		return backb_listOlists

	if file_path[-5:] == ".gff3":
		print('GFF3 detected')
		list_of_dicts = []
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if not line.startswith('#') and len(line) > 1:
				entries = line.split('\t')
				entries[8] = entries[8].strip('\n')
				entries_extra_info = entries[8].split(';')

				if entries[2] == 'gene':

					NOTE = ''
					for add_info in entries_extra_info:
						if 'locus_tag' in add_info:
							LOCUS = add_info[10:]
						if 'Name' in add_info:
							SYMBOL = add_info[5:]
						if 'Note' in add_info:
							NOTE = add_info[5:]

					# row['LOCUS'] row['SYMBOL'] row['SYNOYM'] row['LENGTH']  row['START'] row['STOP']  row['STRAND']  row['NAME']  row['CHROMOSOME']  row['GENOME ONTOLOGY']  row['ENZYME CODE']  row['KEGG']  row['PATHWAY']  row['REACTION'] row['COG'] row['PFAM']  row['OPERON']

					entry_dict = {'CHROM': entries[0], 'LOCUS': LOCUS, 'START': entries[3], 'STOP': entries[4],
								  'STRAND': entries[6], 'SYMBOL': SYMBOL, 'INFO': entries_extra_info, 'NOTE': NOTE}
					list_of_dicts.append(entry_dict)
		return list_of_dicts


def gff2network(inGFF_file, outGML_location):
	"""
	This function takes the uploaded GF file and creates a local genome file in GML format.
	:param inGFF_file: Path to the uploaded GFF file
	:return:
	"""
	#print('Converting')
	#print(inGFF_file)

	GF_dict = input_parser(inGFF_file)

	#print(GF_dict[0])
	#print(GF_dict[-1])

	# Creating the graph
	organism_network = nx.MultiDiGraph(species_ID="species_ID_placeholder", isolate="isolate_placeholder", organism="organism_placeholder")

	# Adding gene nodes
	for entry in GF_dict:
		organism_network.add_node(entry['SYMBOL'], name=entry['SYMBOL'], start=entry['START'], stop=entry['STOP'], strand=entry['STRAND'], type='gene')

	list_pos = 0
	for entry in GF_dict:
		if list_pos != 0:
			organism_network.add_edge(GF_dict[list_pos - 1]['SYMBOL'], entry['SYMBOL'])
		list_pos += 1

	organism_network.add_edge(GF_dict[-1]['SYMBOL'], GF_dict[0]['SYMBOL'])

	# Linking gene nodes

	nx.write_graphml(organism_network, outGML_location + '/test_result_network.xml')

	print("Save the GML file")


def extract_subview_genome(graph_obj, gene_name, link_range=2):

	#node_list = graph_obj.neighbors(gene_name)
	node_list = [gene_name]

	edge_list = graph_obj.out_edges(gene_name) + graph_obj.in_edges(gene_name)
	print(link_range)

	step_number = 1

	while step_number <= link_range:

		for a_gene in node_list:
			print(a_gene)
			for neighbour_node in graph_obj.out_edges(a_gene):
				print(neighbour_node)
				node_list = node_list + [neighbour_node[1]]

			for neighbour_node in graph_obj.in_edges(a_gene):
				print(neighbour_node)
				node_list = node_list + [neighbour_node[0]]

		step_number += 1
		print('ping')

	print(node_list)

	sub_graph = graph_obj.subgraph(node_list)

	return sub_graph


def expression2network(in_gff, in_network, outGML_location):
	#print(os.getcwd())
	expression_list_o_dicts = input_parser(in_gff)
	#print(expression_list_o_dicts[0])

	#print(in_network.node['dnaA'])

	condition_list = []
	individual_condition_list = []

	# Adding the expression data

	for expression_line in expression_list_o_dicts:
		current_gene = expression_line['gene']
		in_network.node[current_gene][expression_line['sample_1'] + '_' + expression_line['sample_2'] + '_' + 'FC'] = expression_line['log2(fold_change)']
		in_network.node[current_gene][expression_line['sample_1'] + '_' + expression_line['sample_2'] + '_' + 'SIG'] = expression_line['significant']
		in_network.node[current_gene][expression_line['sample_1'] + '_' + expression_line['sample_2'] + '_' + 'qVAL'] = expression_line['q_value']
		in_network.node[current_gene][expression_line['sample_1'] + '_' + expression_line['sample_2'] + '_' + 'pVAL'] = expression_line['p_value']

		if (expression_line['sample_1'] + '_' + expression_line['sample_2']) not in condition_list:
			condition_list = condition_list + [expression_line['sample_1'] + '_' + expression_line['sample_2']]

		if expression_line['sample_1'] not in individual_condition_list:
			individual_condition_list.append(expression_line['sample_1'])

		if expression_line['sample_2'] not in individual_condition_list:
			individual_condition_list.append(expression_line['sample_2'])

	in_network.graph['comparisons'] = ",".join(condition_list)
	in_network.graph['conditions'] = ",".join(individual_condition_list)

	#print(in_network.node['dnaA'])

	nx.write_graphml(in_network, outGML_location + '/test_exp_network.xml')

	print("Save the GML file")


def filter_variants(var_object, q_threshold):
	pass_list = []

	for variant in var_object:
		if float(variant['QUAL']) > q_threshold:
			pass_list.append(variant)

	return pass_list


def variants2network(in_vcf, in_network, out_gml_location):

	var_lod = input_parser(in_vcf)

	print(var_lod[0])
	print(len(var_lod))

	filtered_variants = filter_variants(var_lod, 100)

	print(len(filtered_variants))

	for variant in filtered_variants:
		for a_node in in_network.nodes():
			if in_network.node[a_node]['type'] == 'gene':
				print(a_node)
				print(in_network.node[a_node]['start'])

				if int(in_network.node[a_node]['start']) < int(variant['POS']) < int(in_network.node[a_node]['stop']):

					var_node_name = variant['POS'] + '_' + variant['REF'] + '>' + variant['ALT']

					# Create a new node
					if var_node_name not in in_network:
						in_network.add_node(var_node_name, ALT=variant['ALT'], REF=variant['REF'], POS=variant['POS'], type='variant')
						in_network.add_edge(var_node_name, a_node)
					else:
						in_network.add_edge(var_node_name, a_node)

	nx.write_graphml(in_network, out_gml_location + '/test_var_network.xml')

	print('here we go')


def get_undirected_neighbours(graph_obj, node_name):
	node_list = []
	node_list = node_list + [node_name]

	for neighbour_node in graph_obj.out_edges(node_name):
		node_list = node_list + [neighbour_node[1]]

	for neighbour_node in graph_obj.in_edges(node_name):
		node_list = node_list + [neighbour_node[0]]

	return node_list


def qval_product(gene_list, graph_obj, condition):

	q_val_product = 1.0

	for a_gene in gene_list:

		if condition + '_qVAL' in graph_obj.node[a_gene].keys():

			q_val_product = q_val_product * float(graph_obj.node[a_gene][condition + '_qVAL'])

	return q_val_product


def extract_subnets(in_graph_obj, condition='all'):

	list_of_sig_genes = []
	list_of_linked_genes = []

	if condition != 'all':

		# get all gene nodes differentially expressed

		for a_node in in_graph_obj.nodes():

			if condition + '_SIG' in in_graph_obj.node[a_node].keys():
				if in_graph_obj.node[a_node][condition + '_SIG'] == 'yes':

					#print(in_graph_obj.node[a_node])

					list_of_sig_genes.append(a_node)

					# Now start extracting subgraphs.

					connected_nodes = get_undirected_neighbours(in_graph_obj, a_node)
					list_of_linked_genes = list_of_linked_genes + connected_nodes

	sub_graph = in_graph_obj.subgraph(list_of_linked_genes)

	undir_subgraph = sub_graph.to_undirected()

	list_of_subgraph_sets = sorted(nx.connected_components(undir_subgraph), key=len, reverse=True)

	list_of_results = []
	count = 1
	for gene_set in list_of_subgraph_sets:

		# Count number of variant nodes
		var_nodes = 0
		for a_node in gene_set:
			if (in_graph_obj.node[a_node]['type']) == 'variant':
				var_nodes += 1

		set_result = {'label': 'Subnet_' + str(count),
					  	'node_count': len(gene_set),
						'qval_prod': qval_product(gene_set, in_graph_obj, condition),
						'nodes': list(gene_set),
					  	'variants': var_nodes
					  }
		list_of_results.append(set_result)
		count += 1

	#nx.write_graphml(sub_graph, 'processed' + '/sig_subnet.xml')

	return list_of_results


def extract_gene_expression_table(in_graph_obj, a_gene):

	#print(in_graph_obj.graph)
	#print(len(in_graph_obj.graph['conditions'].split(',')))

	node_dict = in_graph_obj.nodes(a_gene)[0][1]

	table_data = []

	for condition in in_graph_obj.graph['comparisons'].split(','):
		cond_entry = {
			'fold_change': float(node_dict[condition + '_FC']),
			'Q_value': float(node_dict[condition + '_qVAL']),
			'Conditions': condition
		}

		table_data.append(cond_entry)

	#print(table_data)

	class GeneTable(tables.Table):
		Conditions = tables.Column()
		fold_change = tables.Column()
		Q_value = tables.Column()

	table = GeneTable(table_data)

	return table


def get_current_settings_dict():

	current_experiment = CurrentSettings.objects.values()[0]['current_experimental_setup']

	experiment_dict = Experiment.objects.all().filter(experiment_shortname=current_experiment)[0].__dict__
	setting_dict = CurrentSettings.objects.values()[0]

	experiment_dict.update(setting_dict)

	return experiment_dict


def get_comparison_list_based_on_settings():

	setting_dict = get_current_settings_dict()
	cond_list = setting_dict['comparisons'].split(',')

	genome_graph = nx.read_graphml('processed/test_var_network.xml')

	comparison_list = genome_graph.graph['comparisons'].split(',')

	filtered_table = []

	for comparison in comparison_list:
		cond_1 = comparison.split('_')[0]
		cond_2 = comparison.split('_')[1]

		cond_1_stripped = re.sub(r'|'.join(map(re.escape, cond_list)), '', cond_1)
		cond_2_stripped = re.sub(r'|'.join(map(re.escape, cond_list)), '', cond_2)

		if cond_1_stripped == cond_2_stripped:
			filtered_table.append(comparison)

	return filtered_table


def extract_and_format_gene_expression_table(in_graph_obj, a_gene):

	#print(in_graph_obj.graph)
	#print(len(in_graph_obj.graph['conditions'].split(',')))

	node_dict = in_graph_obj.node[a_gene]

	table_data = []

	for condition in in_graph_obj.graph['comparisons'].split(','):
		cond_entry = {
			'fold_change': float(node_dict[condition + '_FC']),
			'Q_value': float(node_dict[condition + '_qVAL']),
			'Conditions': condition
		}

		table_data.append(cond_entry)

	current_exp_setup = get_current_settings_dict()

	# Filter based on current experimental design

	cond_list = current_exp_setup['comparisons'].split(',')

	filtered_table = []

	for row in table_data:
		cond_1 = row['Conditions'].split('_')[0]
		cond_2 = row['Conditions'].split('_')[1]

		cond_1_stripped = re.sub(r'|'.join(map(re.escape, cond_list)), '', cond_1)
		cond_2_stripped = re.sub(r'|'.join(map(re.escape, cond_list)), '', cond_2)

		if cond_1_stripped == cond_2_stripped:
			filtered_table.append(row)

	dirty_table = '<table> \n'

	dirty_table = dirty_table + '<tr> <th>Conditions</th> <th>Q value</th> <th>Fold change</th> </tr>'

	for row in filtered_table:
		if row['Q_value'] < 0.05:

			row_string = '<tr bgcolor="#FF0000"> \n' + \
				'<td>' + str(row['Conditions']) + '</td> \n' + \
				'<td>' + str(row['Q_value']) + '</td> \n' + \
				'<td>' + str(row['fold_change']) + '</td> \n' + \
				'</tr>'

			dirty_table = dirty_table + row_string

		else:

			row_string = '<tr> \n' + \
				'<td>' + str(row['Conditions']) + '</td> \n' + \
				'<td>' + str(row['Q_value']) + '</td> \n' + \
				'<td>' + str(row['fold_change']) + '</td> \n' + \
				'</tr>'

			dirty_table = dirty_table + row_string

	dirty_table = dirty_table + '</table>'

	return dirty_table


def subnet_to_json_dict(in_subnet):

	# Get all gene nodes:

	nx.write_graphml(in_subnet, "test.graphml")

	comparison_list = get_comparison_list_based_on_settings()

	result_dict = {}

	node_list = []
	ordered_gene_list = []
	node_data = []
	variant_nodes = []

	# Get a list of genes ordered by chrom position, excluding variants.
	edge_node = []
	for a_node in in_subnet:

		# Find start / stop nodes.
		if in_subnet.degree(a_node) == 1 and in_subnet.node[a_node]['type'] == 'gene':
			edge_node.append(a_node)
		elif in_subnet.node[a_node]['type'] == 'variant':
			variant_nodes.append(a_node)

	ordered_gene_list = nx.shortest_path(in_subnet,source=edge_node[0],target=edge_node[1])

	# Base info
	for comparison in comparison_list:

		base_dict = {
			'label': comparison,
			'backgroundColor': 'grey',
			'borderColor': 'black',
			'borderWidth': 1,
			'data': []
		}

		node_data.append(base_dict)

	for a_node in ordered_gene_list:

		#print(a_node)

		node_list.append(a_node)

		for a_condition in node_data:
			local_cond = a_condition['label']
			if in_subnet.node[a_node]['type'] == 'gene':
				a_condition['data'].append(float(in_subnet.node[a_node][local_cond + '_FC']))
				if in_subnet.node[a_node][local_cond + '_SIG'] == 'yes':
					a_condition['backgroundColor'] = 'red'

	result_dict['labels'] = node_list
	result_dict['datasets'] = node_data

	print('this one')
	print(result_dict)
	print('--------')

	return result_dict


def classify_variant(sequence, seq_start, variant_string, strand):

	# How to deal with 840213_GC>GCTGTTC,GGCTGTTC,GTTC ???

	# 1634523_C>T
	var_position = int(variant_string.split('_')[0])
	var_alt = variant_string.split('>')[1]
	var_ref = variant_string.split('>')[0].split('_')[1]
	print(var_ref)
	slice_start = var_position - seq_start
	var_stop = slice_start + len(var_ref)
	alt_sequence = sequence[:slice_start] + var_alt + sequence[var_stop:]
	print('alt seq')

	if strand == '+':
		aa_sequence = translate(sequence)
		alt_aa_sequence = translate(alt_sequence)

	else:
		rev_seq = reverse_compliment(sequence)
		aa_sequence = translate(rev_seq)

		alt_rev_seq = reverse_compliment(alt_sequence)
		alt_aa_sequence = translate(alt_rev_seq)

	# First, is it synonymous?
	if aa_sequence == alt_aa_sequence:
		mut_classification = 'silent|low'

	elif aa_sequence.count('*') != alt_aa_sequence.count('*'):
		mut_classification = 'New stop codon|high'

	elif len(aa_sequence) == len(alt_aa_sequence):
		mut_classification = "non-synon, inframe|medium"

	elif len(var_ref) % 3 is True and len(var_alt) % 3 is True:
		mut_classification = "non-synon, inframe, insertion|medium"

	else:
		mut_classification = "non-synon, frameshift|medium"

	return mut_classification


def get_variant_info(in_graph_obj, in_gengraph_obj):



	# Get all variants

	var_list = []

	for a_node in in_graph_obj.nodes():

		if in_graph_obj.node[a_node]['type'] == "variant":

			print(a_node)
			var_pos = int(in_graph_obj.node[a_node]['POS'])

			var_list.append(a_node)

			node_neighbours = get_undirected_neighbours(in_graph_obj, a_node)

			# Get genes that contain variants

			for a_neighbour in node_neighbours:
				if in_graph_obj.node[a_neighbour]['type'] == "gene":
					gene_start = int(in_graph_obj.node[a_neighbour]["start"])
					gene_stop = int(in_graph_obj.node[a_neighbour]["stop"])
					if gene_start < var_pos < gene_stop:
						print(a_neighbour)

						gene_seq = extract_original_seq_region_fast(in_gengraph_obj, gene_start, gene_stop, 'H37Rv')

						print(in_graph_obj.node[a_neighbour]["strand"])

						variant_class = classify_variant(gene_seq, gene_start, a_node, in_graph_obj.node[a_neighbour]["strand"])
						print(variant_class)
						print('\n')








	# Determine the effect on the gene

	print('here')

	return 'sdf'


def check_loaded_genomes(dir):

	found_files = {'genome': 'none', 'annotation': 'none'}

	onlyfiles = [f for f in listdir(dir) if isfile(join(dir, f))]

	annotation_filetypes = ['gff', 'gtf', 'bed']

	genome_filetypes = ['fna', 'fa', 'fasta']

	# Find annotations
	for file_thing in onlyfiles:
		extension = file_thing.split('.')[-1]
		if extension in annotation_filetypes:
			found_files['annotation'] = file_thing
		elif extension in genome_filetypes:
			found_files['genome'] = file_thing

	return found_files


def check_loaded_annotations(dir):

	found_files = {'custom_annotation': 'none', 'gene_list': 'none'}

	onlyfiles = [f for f in listdir(dir) if isfile(join(dir, f))]

	annotation_filetypes = ['gff', 'gtf', 'bed', 'gff3']

	# Find annotations
	for file_thing in onlyfiles:
		extension = file_thing.split('.')[-1]
		extension = extension.lower()

		if extension in annotation_filetypes and 'custom_subset.gff' not in onlyfiles:
			found_files['custom_annotation'] = file_thing
		elif extension == 'txt':
			found_files['gene_list'] = file_thing
		elif 'custom_subset.gff' in onlyfiles:
			found_files['custom_annotation'] = 'custom_subset.gff'

	return found_files


def extract_gff_subset(anno_file_path, tar_feat_path, anno_folder):

	tar_feat_file = open(tar_feat_path, 'r')

	search_gene_list = []

	for line in tar_feat_file:
		if len(line) > 1:
			search_gene_list.append(line.strip())

	print(search_gene_list)

	anno_file = open(anno_file_path, 'r')

	out_anno_file = open(anno_folder + "custom_subset.gff", 'w')

	for line in anno_file:
		if line[0] == '#':
			out_anno_file.write(line)

	for a_gene in search_gene_list:

		anno_file = open(anno_file_path, 'r')

		for line in anno_file:

			if line[0] != '#':

				line_list = line.split('\t')

				info_list = line_list[8]

				if a_gene in info_list:

					out_anno_file.write(line)

	anno_file.close()
	out_anno_file.close()

	return 'Done'

