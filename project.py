#Copyright (C) 2017-2018 Alexandros Kanterakis, mail:kantale@ics.forth.gr
###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389

#########################----IMPORTS-------------------------------

import argparse
import gzip
import re
import numpy as np
import pandas as pd
from collections import Counter
import random


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
		return df[df['INFO'].str.contains("VT=SNP")]   ###dataframe contains only SNIP's
	elif re.match('.+vcf$', file_name):
		with open(file_name, 'r') as f:
			lines = [l for l in f if not l.startswith('##')]
		lines=list(map(lambda each:each.replace("\n",""), lines))
		df=pd.Series(lines)
		df=df.str.split("\t", expand=True)
		df.columns = df.iloc[0]
		df=df.reindex(df.index.drop(0))
		return df[df['INFO'].str.contains("VT=SNP")]
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
		z[ z == '0|0' ] = 0
		z[ z == '1|0' ] = 1
		z[ z == '0|1' ] = 1
		z[ z == '1|1' ] = 2
	else:
		print( "Population name: {} was not in the list of names: {}".format(pop.name, pop_dict.keys()))
		assert( 0 )
	return np.array( z.sum( axis = 0 )/z.shape[0] )
	

def create_output_list( file_name, sample_file, pop_name, sample_number ,snps, index) :
	df=read_vcf_file(file_name)
	print ((pop_name))
	frequency=generate_pop_freq(df, sample_file, pop_name)
	header=[pop_name+str(k) for k in range(1,sample_number+1)]
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
	



def simulate(input_file, sample_file, pop_name, snps, my_output_file) :
	'''Both part 3 and 4'''
	sample_list=[]
	data=read_vcf_file(input_file)
	snp_freq_index=[]
	for i in range(0,snps) :
		snp_freq_index.append(random.choice(range(0,data.shape[0])))  ####random snp index 
	for i in range(0,len(pop_name)) :
		sample_list.append(create_output_list(input_file, sample_file, pop_name[i][0], int(pop_name[i][1]), snps, snp_freq_index))

	
	pd.concat(sample_list, axis=1).to_csv(my_output_file, sep="\t", mode='w', index=False)

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






