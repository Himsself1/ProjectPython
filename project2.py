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
import sys
import collections as cl
##################

##################### ------------our arguments-----------------------###############
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', type=str)
parser.add_argument('--action')
parser.add_argument('--sample_filename', type=str)
parser.add_argument('--population', nargs = 2, action='append')
parser.add_argument('--output', type=str)
parser.add_argument('--SNPs', type=int)
parser.add_argument('--independent', type=int)
parser.add_argument('--input_filename', type=str)
parser.add_argument('--PCA_filename', type=str)
parser.add_argument('--PCA_plot', type=str)
parser.add_argument('--ITERATIONS', type=int)
parser.add_argument('---MINIMUM_AF', type=int)
parser.add_argument('--START', type=int)
parser.add_argument('--END', type=int)




args = parser.parse_args()

def read_vcf_file(file_name):
	'''
	Create a pandas table of a .vcf file (even if it's compressed)
	'''
	if file_name.endswith(".vcf.gz") or file_name.endswith(".vcf") :
		comms = 0
		f = gzip.open if file_name.endswith('.gz') else open
		with f(file_name) as file: ####count how many lines are comments
			for line in file:
				line=line.decode() if type(line)==bytes else line
				if line.startswith('##'):
					comms += 1
				else:
					break
		comp = 'gzip' if file_name.endswith('.gz') else None
		df=pd.read_table(file_name, compression=comp, skiprows=comms, header=0)
		df=df[df['INFO'].str.contains("VT=SNP")] 
		return df[~df['INFO'].str.contains("MULTI")]  ###dataframe contains only SNP's
	else :
		raise Exception("Invalid File Extension. File needs to end with .vcf or .vcf.gz")
	
def vcf_info(file_name) :
	df=read_vcf_file(file_name)
	print ("File has", df.shape[1]-9,"samples","\n","File has",df.shape[0],"SNPs" )
	
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

def convert_genotype_to_number(z):
	'''
	Converts a VCF genotype to numerical (zero one or two mutations per sample per SNP position) 
	'''
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 2
	return z

def generate_pop_freq( my_data, sample_data, pop_name, independent=None ):
	'''
	Generates a frequency vector for each SNP for individuals in a given population
	'''
	## Makes a dictionary with keys = population names and values = individuals in each population
	pop_dict = dict()
	for i in range(0, sample_data.shape[0]):
		if sample_data[i,1] in pop_dict.keys():
			pop_dict[sample_data[i,1]] = pop_dict[sample_data[i,1]] + [sample_data[i,0]]
		else:
			pop_dict[sample_data[i,1]] = [sample_data[i,0]]
	## Swap genotypes with numbers to make calculations easier
	if pop_name in pop_dict.keys() :
		z = np.array([my_data[x] for x in pop_dict[pop_name]])
		ret=convert_genotype_to_number(z)
	else:
		raise Exception( "Population name: {} was not in the list of names: {}".format(pop.name, pop_dict.keys()))
	return np.array( ret.sum( axis = 0 )/ (2*ret.shape[0]) )

def create_output_list( frequency, pop_name, sample_number , index) :
	'''
	Creates the simulated data for one population
	'''
	output_list = []
	for it in range(0,len(index)):
		snp_freq = frequency[index[it]]
		ones = int(snp_freq*2*float(sample_number)) ## Posa 1 8a uparxoun sto dataset? Takes the integer part of the number
		if ones==0 :
			genotype = ["0"+'|'+"0" for k in range(0,sample_number*2-1,2) ] ## Kanw concatenate to genotype
			output_list.append(genotype)
		else :
			a = random.sample(range(2*sample_number), ones)
			output = [0]*2*sample_number
			for i in a:
				output[i] = 1 ## Ftiaxnw tous assous
			genotype = [str(output[k])+'|'+str(output[k+1]) for k in range(0,len(output)-1,2) ]
			output_list.append(genotype)
	
	return pd.DataFrame(output_list)

def simulation(data, sample_table, population_list, snps,  independent=None) :
	populations=[i[0] for i in population_list]
	population_sizes={i[0]:int(i[1]) for i in population_list}
	header=[pop_name+"_"+str(k) for pop_name in populations for k in range(1,int(population_sizes[pop_name])+1)  ]
	population_frequences={i:generate_pop_freq(data,sample_table, i ) for i in populations}
	snp_freq_index=[random.choice(range(0,data.shape[0])) for i in range(0,snps)]
	output_list=[create_output_list(population_frequences[k], k , population_sizes[k], snp_freq_index) for k in populations]
 
	if snps > 0 and independent==None:
		dependent_frame=pd.concat(output_list, axis=1)
		dependent_frame.columns=header
		return dependent_frame
	elif snps > 0 and independent>0 :
		dependent_frame=pd.concat(output_list, axis=1)
		dependent_frame.columns=header
		independent_frequencies = np.array([ float(population_sizes[i])*population_frequences[i] for i in populations ]).sum(axis = 0)/float(sum(population_sizes.values()))
		independent_output_list=create_output_list(independent_frequencies, "independent" , sum(population_sizes.values()),[random.choice(range(0,data.shape[0])) for i in range(0,independent)])
		independent_output_list.columns=header
		return pd.concat((dependent_frame,independent_output_list), axis=0) 
	elif snps == 0 and  independent>0:
		independent_frequencies = np.array([ float(population_sizes[i])*population_frequences[i] for i in populations ]).sum(axis = 0)/float(sum(population_sizes.values()))
		independent_output_list=create_output_list(independent_frequencies, "independent" , sum(population_sizes.values()),[random.choice(range(0,data.shape[0])) for i in range(0,independent)])
		independent_output_list.columns=header
		return independent_output_list
	else :
		raise Exception ('Independent variable and SNPs number must be zero or positive integers only' )

def read_my_labels( my_vector ):
	'''
	Takes a pd.Series with label names and returns grouped labels
	'''
	kefali = pd.DataFrame(my_vector.str.split("_").tolist(), columns=["population", "index"])
	kefali["index"] = pd.to_numeric(kefali['index'])
	aa = kefali.groupby("population").count()
	return aa

def my_pca(np_table) :
	'''
	Principal component analysis to simulated data tables
	'''
	genotypes_array=convert_genotype_to_number(np_table)
	pca = PCA(n_components=2)
	pca.fit(genotypes_array.T)
	genotypes_PCA = pca.transform(genotypes_array.T)
	return genotypes_PCA
	
def pca_plot(genotypes_PCA, header) :
	'''
	Saves two-dimensional data produced b PCA
	'''	
	aa = read_my_labels(header)
	left = 0
	for i in aa.values:
		plt.plot(genotypes_PCA[left:left+i[0],0], genotypes_PCA[left:left+i[0],1], '.', label=i)
		left = left + i[0]
	plt.savefig(args.PCA_plot)

def find_majority( mia_lista, error, exclude ):
	'''
	Briskei to stoixeio me to megalutero frequency sth mia_lista
	mono ean den to exei 3anabrei
	'''
	count_dict=cl.Counter( mia_lista )
	b=(max( [[count_dict[kleidi], kleidi] for kleidi in count_dict.keys() if kleidi not in exclude]))
	error.append( mia_lista.shape[0] - mia_lista[ mia_lista == b[1] ].shape[0] )
	exclude.append(b[1])
	return error, exclude

def k_means(np_table):
	'''
	Takes as input the output of my_pca function!
	Kanei k-means analoga me to posoi diaforetikoi plh8usmoi uparxoun sto deigma
	Kai briskei to error rate
	'''
	header=pd.Series(np_table[:,0])
	aa = read_my_labels(header)
	genotypes_PCA=np_table[:,1:].astype(float)
	kmeans = KMeans(n_clusters=aa['index'].shape[0], random_state=0).fit(genotypes_PCA)
	left = 0
	error = []
	exclude = []
	for i in aa['index']:
		error, exclude = find_majority( kmeans.labels_[left: left+i], error, exclude )
		left += i
	total_error = sum(error)/kmeans.labels_.shape[0]
	return total_error
 
def find_ratio(data, sample_table, population_list, independent=None):
	simulation_table=simulation(data, sample_table, population_list, 0, independent)
	simulation_table=np.vstack((np.array(simulation_table.columns),np.array(simulation_table)))
	pca_table=my_pca(simulation_table[1:,])
	k_means_value = k_means(np.column_stack((simulation_table[0,:],pca_table)))
	while True :                                                    ################<<<<<<---- if's need work
		dependent_simulation_table=simulation(data, sample_table, population_list, 10)
		simulation_table=np.vstack((simulation_table,dependent_simulation_table))
		pca_table=my_pca(simulation_table[1:,])
		if abs(k_means_value - k_means(np.column_stack((simulation_table[0,:],pca_table))))>0.01 and k_means(np.column_stack((simulation_table[0,:],pca_table)))>0.2 :
			k_means_value=k_means(np.column_stack((simulation_table[0,:],pca_table)))
		elif k_means(np.column_stack((simulation_table[0,:],pca_table)))<0.2 :
			#return k_means(np.column_stack((simulation_table[0,:],pca_table)))
			return simulation_table.shape[0]-independent-1
			break

	



################------- THIS IS HOW WE DO IT -----------################

action=args.action
if action == "VCF_INFO" :   ###part1
	vcf_info(args.vcf)
elif action == "VALIDATE_SAMPLE_INFO" : ####part 2
	validate_sample(args.vcf, args.sample_filename)
elif action== "SAMPLE_INFO" :   ####part 2
	sample_info(args.sample_filename)
elif action == "SIMULATE" :  ####parts 3-5
	data=read_vcf_file(args.vcf)
	sample_table = np.loadtxt( args.sample_filename, delimiter = "\t", skiprows = 1, dtype = "object" )  ## Stores the decription of individuals
	simulation_table=simulation(data, sample_table , args.population, args.SNPs,  args.independent)
	simulation_table.to_csv(args.output, sep="\t", mode='w', index=False)
elif action == "PCA" : ### part 6
	np_table=np.loadtxt(args.input_filename, delimiter="\t", dtype="object")
	header=pd.Series(np_table[0,])
	pca_table=my_pca(np_table[1:,])
	pca_plot(pca_table, header)  ### for the plot 
	header=np.array(header)
	np.savetxt(args.PCA_filename, np.column_stack((header,pca_table)) ,fmt="%s\t%f\t%f")
elif action == "CLUSTER" : ##part 7
	pca_table=np.loadtxt(args.PCA_filename, delimiter = '\t', dtype = 'object')
	print(k_means(pca_table))
elif action == "FIND_RATIO" : ##part8
	data=read_vcf_file(args.vcf)
	sample_table = np.loadtxt( args.sample_filename, delimiter = "\t", skiprows = 1, dtype = "object" )  ## Stores the decription of individuals
	ratio=find_ratio(data, sample_table, args.population ,args.independent)
	print (ratio)
	
