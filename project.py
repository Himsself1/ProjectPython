#Copyright (C) 2017-2018 Alexandros Kanterakis, mail:kantale@ics.forth.gr
###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389

#########################----IMPORTS-------------------------------

import argparse
import gzip
import re
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import random
import matplotlib.pyplot as plt
import matplotlib as mpl



##################### ------------our arguments-----------------------
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', nargs = 1)  ###θα παίρνει ένα vcf ή ένα vcf.gz αρχείο
parser.add_argument('--action', nargs = 1)  ##Το πρόγραμμά σας θα τυπώνει πόσα SNPs έχει το vcf και πόσα samples.



input_file = 'chr22_25_lines.vcf'
sample_file= 'sample_information.csv'
## input_file = args.vcf

###########------------our functions------------------

def read_vcf_file(file_name):
	'''
	Create a pandas table with header 
	'''
	if re.match('.+vcf.gz$', file_name):
		with gzopen(file_name, 'r') as f:
			lines = [l for l in f if not l.startswith('##')]
		lines=list(map(lambda each:each.replace("\n",""), lines))
		df=pd.Series(lines)
		df=df.str.split("\t", expand=True)
		df.columns = df.iloc[0]
		df=df.reindex(df.index.drop(0))
		df=df[df['INFO'].str.contains("VT=SNP")] 
		return df[~df['INFO'].str.contains("MULTI")]  ###dataframe contains only SNP's
	elif re.match('.+vcf$', file_name):
		with open(file_name, 'r') as f:
			lines = [l for l in f if not l.startswith('##')]
		lines=list(map(lambda each:each.replace("\n",""), lines))
		df=pd.Series(lines)
		df=df.str.split("\t", expand=True)
		df.columns = df.iloc[0]
		df=df.reindex(df.index.drop(0))
		df=df[df['INFO'].str.contains("VT=SNP")] 
		return df[~df['INFO'].str.contains("MULTI")] 
		
	else :
		raise Exception("Invalid File Extension")
	
def vcf_info(file_name) :
	df=read_vcf_file(file_name)
	print ("File has", df.shape[1]-9,"samples" ,"\n","File has",df.shape[0],"SNPs" )


def sample_info(sample_file):
	df=pd.read_csv(sample_file, sep="\t", header = 0)
	print ("File has", len(df.super_pop.unique()), "Areas.")
	for i in range(0,len(df.super_pop.unique())) :
		print("Area", i+1, "is",
			df.super_pop.unique()[i],
			"and contains",
			df[(df["super_pop"]==df.super_pop.unique()[i])].shape[0],
			"samples splitted in the following populations:")
		q=df[(df["super_pop"]==df.super_pop.unique()[i])]
		for k in range(0,len(q['pop'].unique())) :
			print (q['pop'].unique()[k], q[(q["pop"]==q['pop'].unique()[k])].shape[0],"samples")
	
def validate_sample(file_name, sample_file ) :
	our_samples=np.array(read_vcf_file(file_name).columns.values[9:])
	id_array=np.array(pd.read_csv(sample_file, sep="\t", header = 0)["sample"].tolist())
	if np.array_equiv(np.sort(our_samples),np.sort(id_array)):
		print ("Everything is OK!")
	else :
		b = np.setdiff1d(id_array,our_samples)
		q = np.setdiff1d(our_samples,id_array)
		if len(b) != 0:
			print( "These samples are present in VCF but not in SAMPLE:", b )
		if len(q) != 0:
			print( "These samples are present in SAMPLE but not in VCF:", q )
			
			
def generate_pop_freq( my_data, sample_file, pop_name ):
	'''
	Generates a frequency vector for each SNP for individuals in a given population
	'''
	## Stores the decription of individuals
	pop_parameters = np.loadtxt( sample_file, delimiter = "\t", skiprows = 1, dtype = "object" ) 
	## Makes a dictionary with keys = population names and values = individuals in each population
	pop_dict = dict()
	for i in range(0, pop_parameters.shape[0]):
		if pop_parameters[i,1] in pop_dict.keys():
			pop_dict[pop_parameters[i,1]] = pop_dict[pop_parameters[i,1]] + [pop_parameters[i,0]]
		else:
			pop_dict[pop_parameters[i,1]] = [pop_parameters[i,0]]
	## Swap genotypes with numbers to m ake calculations easier
	if pop_name in pop_dict.keys():
		z = np.array([my_data[x] for x in pop_dict[pop_name]])
		ret=convert_genotype_to_number(z)
	
	else:
		print( "Population name: {} was not in the list of names: {}".format(pop.name, pop_dict.keys()))
		assert( 0 )
	## print (ret[0:10,])
	return np.array( ret.sum( axis = 0 )/ret.shape[0] )

def convert_genotype_to_number(z):
	'''
	Oti leei to onoma.
	'''
	if z[ z not in np.array(['0|0','1|0','0|1','1|1',])].shape[0] != 0:
		## print( z[ z not in np.array(['0|0','1|0','0|1','1|1',])] )
		## raise Exception("INVALID_SNP_GENOTYPE")
		pass
	
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 2
	return z
	
def generate_independent_freq( my_data, sample_file, pop_nams ):
	'''
	Generates a frequency vector for each SNP for individuals in a given population
	'''
	## Stores the decription of individuals
	pop_parameters = np.loadtxt( sample_file, delimiter = "\t", skiprows = 1, dtype = "object" ) 
	## Makes a dictionary with keys = population names and values = individuals in each population
	pop_dict = dict()
	for i in range(0, pop_parameters.shape[0]):
		if pop_parameters[i,1] in pop_dict.keys():
			pop_dict[pop_parameters[i,1]] = pop_dict[pop_parameters[i,1]] + [pop_parameters[i,0]]
		else:
			pop_dict[pop_parameters[i,1]] = [pop_parameters[i,0]]
	## Swap genotypes with numbers to m ake calculations easier
	multiple_pops=[]
	for pop_name in pop_nams :
		if pop_name in pop_dict.keys():
			multiple_pops.append(np.array([my_data[x] for x in pop_dict[pop_name]]))
		else:
			print( "Population name: {} was not in the list of names: {}".format(pop.name, pop_dict.keys()))
			assert( 0 )
	z=np.concatenate(multiple_pops)
	ret=convert_genotype_to_number(z)
	return np.array( ret.sum( axis = 0 )/ret.shape[0] )
	
	
def create_output_list( file_name, sample_file, pop_name, sample_number ,snps, index) :
	df=read_vcf_file(file_name)
	frequency=generate_pop_freq(df, sample_file, pop_name)
	header=[pop_name+"_"+str(k) for k in range(1,sample_number+1)]
	output_list = []
	for it in range(0,snps):
		snp_freq = frequency[index[it]]
		ones = int(snp_freq*2*sample_number // 1) ## Posa 1 8a uparxoun sto dataset? Takes the integer part of the number
		if ones==0 :
			genotype = []
			for k in range(0,sample_number*2-1,2):
				genotype.append("0"+'|'+"0") ## Kanw concatenate to genotype
			output_list.append(genotype)
		else :
			a = random.sample(range(2*sample_number), ones)
			output = [0]*2*sample_number
			for i in a:
				output[i] = 1 ## Ftiaxnw tous assous
			genotype = []
			for k in range(0,len(output)-1,2):
				genotype.append(str(output[k])+'|'+str(output[k+1])) ## Kanw concatenate to genotype
			output_list.append(genotype)
	return pd.DataFrame(output_list, columns=header)
	
def create_output_independent_list( file_name, sample_file, pop_names, sample_number ,snps, index) :
	df=read_vcf_file(file_name)
	frequency=generate_independent_freq(df, sample_file, pop_names)
	header=["Independent_"+str(k) for k in range(1,sample_number+1)]
	output_list = []
	for it in range(0,snps):
		snp_freq = frequency[index[it]]
		ones = int(snp_freq*2*sample_number // 1) ## Posa 1 8a uparxoun sto dataset? Takes the integer part of the number
		if ones==0 :
			genotype = []
			for k in range(0,sample_number*2-1,2):
				genotype.append("0"+'|'+"0") ## Kanw concatenate to genotype
			output_list.append(genotype)
		else :
			a = random.sample(range(2*sample_number), ones)
			output = [0]*2*sample_number
			for i in a:
				output[i] = 1 ## Ftiaxnw tous assous
			genotype = []
			for k in range(0,len(output)-1,2):
				genotype.append(str(output[k])+'|'+str(output[k+1])) ## Kanw concatenate to genotype
			output_list.append(genotype)
	return pd.DataFrame(output_list, columns=header)


def simulate(input_file, sample_file, pop_name, snps, my_output_file, independent=None) :
	'''Both part 3 and 4'''
	sample_list=[]
	data=read_vcf_file(input_file)
	snp_freq_index=[]
	pop_names=[]
	for i in range(0,snps) :
		snp_freq_index.append(random.choice(range(0,data.shape[0])))  ####random snp index 
	for i in range(0,len(pop_name)) :
		pop_names.append(pop_name[i][0])
		sample_list.append(create_output_list(input_file, sample_file, pop_name[i][0], int(pop_name[i][1]), snps, snp_freq_index))
	if independent==None :
		pd.concat(sample_list, axis=1).to_csv(my_output_file, sep="\t", mode='w', index=False)
	else :
		sample_list.append(create_output_independent_list(input_file, sample_file, pop_names, independent, snps, snp_freq_index))
		
	pd.concat(sample_list, axis=1).to_csv(my_output_file, sep="\t", mode='w', index=False)
	
	
def my_pca(input_file, output_file, plot_file) :
	np_table=np.loadtxt(input_file, delimiter="\t", dtype="object")
	header=pd.Series(np_table[0,])
	aa = read_my_labels(header)
	genotypes_array=convert_genotype_to_number(np_table[1:,])
	pca = PCA(n_components=2)
	pca.fit(genotypes_array.T)
	genotypes_PCA = pca.transform(genotypes_array.T)
	left = 0
	for i in aa.values:
		plt.plot(genotypes_PCA[left:left+i[0],0], genotypes_PCA[left:left+i[0],1], '.', label=i)
		left = left + i[0]
	## plt.plot(genotypes_PCA[100:,0], genotypes_PCA[100:,1], '.', color="blue")
	plt.savefig(plot_file)
	header=np.array(header)
	fmt="%s\t%f\t%f"
	q=np.column_stack((header,genotypes_PCA))
	np.savetxt(output_file, q ,fmt=fmt)
#####---------PART1-----------

#vcf_info(input_file)

#####---------PART2-------------


#sample_info(sample_file)
#validate_sample(input_file, sample_file)

#####--------PART3 and 4-----------

pop_name="GBR"
sample_number=10
snps=100



parser.add_argument('--population', nargs = 2, action='append')  ###θα παίρνει ένα vcf ή ένα vcf.gz αρχείο

args = parser.parse_args()
pop_nams=args.population
#simulate(input_file, sample_file, pop_nams, snps, "output.vcf")


####---------PART5----------------------
#simulate(input_file, sample_file, pop_nams,snps, "output.vcf", 10)

####----------PART6------------------------------

###PCA 


	
	
##my_pca("output.vcf", "pca_file.tsv", "pca_plot.pdf")	


######---------PART7-------------------------

def k_means(input_file) :
        '''
        Takes as input the output of my_pca function!
        '''
	np_table=np.loadtxt(input_file, delimiter="\t", dtype="object")
	header=pd.Series(np_table[:,0])
	aa = read_my_labels(header)
	genotypes_PCA=np_table[:,1:]
	kmeans = KMeans(n_clusters=4, random_state=0).fit(genotypes_PCA)
	print (kmeans.labels_)
	# belong_to_1_clustered_2 = sum([1 for i,x in enumerate(kmeans.labels_) if i<aa["YRI"] and x==0])
	# belong_to_2_clustered_1 = sum([1 for i,x in enumerate(kmeans.labels_) if i>=aa["YRI"] and x==1])
	# errors = belong_to_1_clustered_2 + belong_to_2_clustered_1
	# total = genotypes_PCA.shape[0]
	# success = (total-errors)/total
	# print (success)
	
def read_my_labels( my_vector ):
	'''
	Takes a pd.Series with label names and returns grouped labels
	'''
	kefali = pd.DataFrame(my_vector.str.split("_").tolist(), columns=["population", "index"])
	kefali["index"] = pd.to_numeric(kefali['index'])
	aa = kefali.groupby("population").count()
	return aa

k_means("pca_file.tsv")



############### Part 8 ###########


