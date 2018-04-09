
import pandas as pd
import numpy as np
import allel

# filename = '1kgp_chr1.vcf'
# file = open(filename, 'rt')
# line = file.readline()
# while line[:2] == '##':
	# line = file.readline()
# header=pd.Series(line.split("\t"))

# file.close()

def population_columns( sample_matrix, pop_name ):
	'''
	Returns the names of all individuals of population <pop_name>
	<sample_matrix> must be a pandas matrix of the sample_information file
	<pop_name> is ONE POPULATION as string
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	if sample_matrix['pop'].str.contains(pop_name).sum() == 0:
		raise Exception ('Error, Wrong population name: {}'.format(pop_name))
	else:
		return sample_matrix[sample_matrix['pop'].str.contains(pop_name)]['sample']
		
def calculate_frequencies( data, individuals ):
	'''
	Returns the frequencies of all snps contained in <data[individuals]>
	<data> is a pandas Data Frame
	<individuals> is a 1D pandas Data Frame with the names of individuals of the <data> 
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	callset = allel.read_vcf(data, samples=individuals ,fields=['calldata/GT'])
	gt=allel.GenotypeDaskArray(callset['calldata/GT'])
	no_alleles=gt.count_alleles().compute()
	frequencies=np.sum(np.array(no_alleles[:,1:]),axis=1)
	snp_type = allel.vcf_to_dataframe(data, fields=[ 'variants/is_snp'])
	snp_only = snp_type[snp_type['is_snp'] == True]
	index = list(snp_only.index.values)
	return pd.DataFrame(frequencies[index])/(2*len(list(individuals)))

sample_table=pd.read_csv("sample_information.csv", sep="\t", header = 0)
population_sizes=dict(sample_table["pop"].value_counts())
populations=list(population_sizes.keys())

pop_frequencies={}
pop_samples={i:population_columns(sample_table,i) for i in populations}
df_list=[]
for q in populations:
	s=calculate_frequencies("chr22_200000_lines.vcf", pop_samples[q])  ###<<<<<<<<<<
	df_list.append(s)
	print(len(df_list),q)
q=pd.concat(df_list, axis=1)
q.columns=populations


q.to_csv("frequencies.csv", sep="\t", index=False)


