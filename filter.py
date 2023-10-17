import os
import pandas as pd


def vcf_filter():
	"""Function filters out list of muttaions that are synonymous and did not pass the VCF"""
	for file in os.listdir():
		if file.endswith('vcf'):
			my_file = open(file).read().split('\n')
			filter_file = open(file.split('.vcf')[0] + '_missense_variants.vcf', 'w')
			for line in my_file:
				if "missense_variant" in line:
					filter_file.write(line + '\n')
			filter_file.close()

def prot_mutation():
	"""Function isolates proteins and muttaions and appends them to a CSV file"""
	for file in os.listdir():
		df =  pd.DataFrame(columns=['Protein','Mutations'])
		if file.endswith('_missense_variants.vcf'):
			my_file = open(file).read().split('\n')
			print(file)
			for line in my_file:
				if line != '':
					if 'p.' in line:
						vcf_line_data = line.split('\t')
						protein = vcf_line_data[7].split('|')[3]
						mutation = vcf_line_data[7].split('|')[-6]
						df = df._append({'Protein':protein, 'Mutations':mutation.split('p.')[1]},ignore_index=True)
			df.to_csv(file.split('.vcf')[0] + '.csv')
			print('\n')

def prot_commons():
	"""Function identifies individual muttaions filtering duplicates amongst samples"""
	prot_mutations_list = []
	prot_mutations_dict = {}
	for file in os.listdir():
		if file.endswith('.csv'):
			data = pd.read_csv(file)
			data.drop('Unnamed: 0', inplace=True, axis=1)
			for ind in data.index:
				if (data['Protein'][ind] + '_' + data['Mutations'][ind]) not in prot_mutations_list:
					prot_mutations_list.append(data['Protein'][ind] + '_' + data['Mutations'][ind])
					prot_mutations_dict.setdefault(data['Protein'][ind], []).append(data['Mutations'][ind])
	prot_df = pd.DataFrame.from_dict(prot_mutations_dict, orient ='index')
	print(prot_df)
	print(len(prot_mutations_list))
	file_name = os.path.basename(os.getcwd())
	prot_df.to_csv(file_name + '.csv')

#vcf_filter()
#prot_mutation()
prot_commons()

