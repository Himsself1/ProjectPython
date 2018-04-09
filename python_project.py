#Copyright (C) 2017-2018 Alexandros Kanterakis, mail:kantale@ics.forth.gr
###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389

#########################----IMPORTS-------------------------------

import argparse
import gzip
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import random
import matplotlib.pyplot as plt
import collections as cl
from sklearn.metrics.pairwise import pairwise_distances
from scipy.cluster import hierarchy
import scipy.stats as st
import math as math
import itertools
import allel
from pathlib import Path

##################

##################### ------------our arguments-----------------------###############
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', type=str, help="Type your .vcf or .vcf.gz file. If your file is in different destination,you need to type the path too")
parser.add_argument('--action', help="Type what action you want this program to perform", choices=["VCF_INFO","SAMPLE_INFO","VALIDATE_SAMPLE_INFO","SIMULATE","PCA","CLUSTER","FIND_RATIO","DENDROGRAM"])
parser.add_argument('--sample_filename', type=str, help="Type the sample file")
parser.add_argument('--population', nargs = 2, action='append', metavar=('name','size'),help="Type the name of the population and the desired sample number")
parser.add_argument('--output', type=str, help="Type the name of the file you want to save the result with the desired file type and destination")
parser.add_argument('--SNPs', type=int, help="Type how many number of SNPs you want to simulate")
parser.add_argument('--independent', type=int,help="Type how many independent SNPs you want to simulate" )
parser.add_argument('--input_filename', type=str, help="Type the name of the file you produced with the simulation action")
parser.add_argument('--PCA_filename', type=str, help="Type the name of the file you want to save your PCA results")
parser.add_argument('--PCA_plot', type=str, help="Type the name of the file you want to save your PCA plot")
parser.add_argument('--iterations', type=int, default = 1, help="Type how many times you want to find the independent-dependent ratio" )
parser.add_argument('--MINIMUM_AF', type=float,default=0, help="Type the minimum frequency of an SNP for the simulation" )
parser.add_argument('--START', type=int, default = 0, help="Type the minimum position of SNPs you want to simulate" )
parser.add_argument('--END', type=int, help="Type the maximum position of SNPs you want to simulate")
parser.add_argument('--jaccard', type=str, help="Type anything to perform an action with Jaccard Index matrices")


args = parser.parse_args()


def vcf_info(file_name, start = 0, end = None) :
	'''
	Prints how many samples and SNPs are in your .vcf file
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	callset=allel.read_vcf(file_name, fields=['samples', 'variants/is_snp'])
	snps_number=np.count_nonzero(callset['variants/is_snp'])
	sample_number=len(callset['samples'])
	print ("File has", sample_number,"samples","\n","File has",snps_number,"SNPs" )

def sample_info(sample_file):
	'''
	Prints how many populations and super-populations (geographic regions) are in the sample information file
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
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

def validate_sample(file_name, sample_file ):
	'''
	Confirms if the sample infomation file and the .vcf file have the same individuals
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	our_samples=np.array(allel.read_vcf(file_name, fields=['samples'])["samples"])
	id_array=np.array(pd.read_csv(sample_file, sep="\t", header = 0)["sample"].tolist())
	if np.array_equiv(np.sort(our_samples),np.sort(id_array)):
		print ("Everything is OK!")
	else:
		b = np.setdiff1d(id_array,our_samples)
		q = np.setdiff1d(our_samples,id_array)
		if len(b) != 0:
			print( "These samples are present in VCF but not in SAMPLE:", b )
		if len(q) != 0:
			print( "These samples are present in SAMPLE but not in VCF:", q )

def convert_genotype_to_number(z):
	'''
	Converts a VCF genotype to numerical (zero one or two mutations per sample per SNP position)
	<z> must be a numpy table or pandas table of genotypes 
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	z = np.array(z)
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 2
	return z

def convert_to_jaccard(z):
	'''
	Converts a VCF genotype to jaccard (zero or one mutations per sample per genotype position)
	<z> must be a numpy table or pandas table of genotypes 
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	z = np.array(z)
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 1
	return z

def jaccard_index (table):
	'''
	calculates the jaccard matrix accordring to hamming distance
	<table> must be a numpy table
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	jac_table = convert_to_jaccard(table)
	return 1 - pairwise_distances(table.T, metric = "hamming")

def find_non_zero_AFs( dict_of_pop_freqs, min_af = 0.0 ):
	'''
	Returns the indexes of acceptable SNPs
	<min_af> is the lowest acceptable frequency
	<dict_of_pop_freqs> is a dictionary with keys = population names and values = frequency vectors for SNPs
	values are numpy arrays
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	table = np.array([ i for i in dict_of_pop_freqs.values() ]).sum( axis = 0 )/len( dict_of_pop_freqs.keys() )
	return np.array( np.where( table >= min_af )[0] ) 

def create_genotype( frequenc, sample_number ):
	'''
	Generates a genotype of size <sample_number> with snp frequency = <frequenc>
	<frequenc> must be a float
	<sample_number> must me float converible
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	ones = int(2*float(sample_number)*frequenc)
	zeros = 2*sample_number - ones
	random_genotypes = [str(0)]*zeros + [str(1)]*ones
	random.shuffle(random_genotypes)
	output_genotype = [(random_genotypes[i] + '|' + random_genotypes[i+1]) for i in range(0, sample_number*2, 2)]
	return output_genotype

def calculate_frequencies( data, individuals ):
	'''
	Returns the frequencies of all snps contained in <data[individuals]>
	<data> is the vcf file containing all the genotypes
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

def create_frequency_file(data, sample_table) :
	'''
	If our file containing all the SNPs frequencies doesnt exists, it creates it
	<data> is the vcf file containing all the genotypes
	<sample_filename> is the sample_information file
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	population_sizes=dict(sample_table["pop"].value_counts())
	populations=list(population_sizes.keys())
	pop_frequencies={}
	pop_samples={i:population_columns(sample_table,i) for i in populations}
	df_list=[]
	for q in populations:
		s=calculate_frequencies("chr22_200000_lines.vcf", pop_samples[q])  ###<<<<<<<<<<
		df_list.append(s)
		print(len(df_list),"Now creating :",q)
	q=pd.concat(df_list, axis=1)
	q.columns=populations
	q.to_csv("frequencies.csv", sep="\t", index=False)

def read_or_create_frequency_file(sample_filename, vcf):
	'''
	Checks if file containing SNPs frequencies exists.  If not it creates it by calling "create_frequency_file" function
	<sample_filename> is the sample_information file
	<vcf> is the file containing all the genotypes
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	sample_table = pd.read_csv(sample_filename, sep="\t", header = 0)
	if Path("frequencies_pop.csv").is_file():
		pass
	else:
		create_frequency_file(vcf, sample_table )
	freq_table=pd.read_csv("frequencies_pop.csv", sep="\t", header=0)
	return freq_table

def simulation_dependent( pop_vecs, pop_names, sample_sizes, total_snps, min_af = 0 ):
	'''
	Simulates dependent genotypes
	<pop_vectors> is the pandas DataFrame containing genotypes frequencies
	<sample_matrix> is a pandas DataFrame containing information from the sample information file
	<pop_names> is a python list with all the names of the populations provided by the command line
	<sample_sizes> is a python dictionary with keys = pop_names and values = how many individuals are to be simulated
	for each population
	<total_snps> is an integer containing the number of snps to be simulated, provided by command line
	<min_af> is the lnowest acceptable frequency of a SNP
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	index=find_non_zero_AFs(pop_vecs, min_af)
	q=np.random.choice(index, size = int(total_snps), replace = True)
	output_dict = { pop_names[i]: pd.DataFrame([create_genotype(pop_vecs[pop_names[i]][ind], sample_sizes[pop_names[i]]) for ind in q],
									columns=[str(pop_names[i])+"_"+str(k) 
									for k in range(1,int(sample_sizes[pop_names[i]])+1)]) 
					for i in range(0,len(pop_names)) }
	finale=pd.concat(output_dict.values(), axis=1)
	return finale

def simulation_independent( pop_vecs, sample_table ,pop_names, sample_sizes, independent, min_af = 0 ):
	'''
	Simulates independent genotypes based on genotype frequencies of all populations in <pop_names>
	<pop_vectors> is the pandas DataFrame containing genotypes frequencies
	<sample_matrix> is a pandas DataFrame containing information from the sample information file
	<pop_names> is a python list with all the names of the populations provided by the command line
	<sample_sizes> is a python dictionary with keys = pop_names and values = how many individuals are to be simulated for each population
	<min_af> is the lnowest acceptable frequency of a SNP
	<independent> is the number of SNPs to be simulated. It is provided by command line
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	population_sizes=dict(sample_table["pop"].value_counts())
	new_pop_vecs={i:pop_vecs[i]*population_sizes[i] for i in pop_vecs.keys()}
	total_freq=sum(new_pop_vecs.values())/sum(population_sizes[k] for k in pop_vecs.keys())
	index=np.array( np.where( total_freq >= min_af )[0] )
	q=np.random.choice(index, size = independent, replace = True)
	output_dict = { pop_names[i]: pd.DataFrame([create_genotype(total_freq[ind], sample_sizes[pop_names[i]]) for ind in q],
									columns=[str(pop_names[i])+"_"+str(k) 
									for k in range(1,int(sample_sizes[pop_names[i]])+1)]) 
				for i in range(0,len(pop_names)) }
	finale=pd.concat(output_dict.values(), axis=1)
	return finale	

def read_my_labels( my_vector ):
	'''
	Takes a pd.Series with label names and returns a dictionary with keys = population labels and values = the indexes of individuals accroding to labeling
	<my_vector> is a pd.Series containing the names of the samples in following format: "POPNAME_INDEX"
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	kefali = pd.DataFrame(my_vector.str.split("_").tolist(), columns=["population", "index"])
	kefali["index"] = pd.to_numeric(kefali['index'])
	aa = { i:[] for i in kefali['population'].unique() }
	[ aa[kefali.as_matrix()[i,0]].append(i) for i in range(0,kefali.shape[0]) ]
	return aa

def my_pca(np_table, jaccard = None) :
	'''
	Principal component analysis to simulated data tables
	<np_table> is a numpy table
	if <jaccard> = True, it performs PCA according to jaccard index
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	if jaccard== None :
		genotypes_array=convert_genotype_to_number(np_table)
		pca = PCA(n_components=2)
		pca.fit(genotypes_array.T)
		genotypes_PCA = pca.transform(genotypes_array.T)
	else:
		genotypes_array=convert_to_jaccard(np_table)
		pca = PCA(n_components=2)
		pca.fit(jaccard_index(np_table))
		genotypes_PCA = pca.transform(jaccard_index(np_table))
	return genotypes_PCA

def pca_plot( genotypes_PCA, header, plot_filename ):
	'''
	Plots two-dimensional data produced by PCA
	<genotypes_PCA> is the output of the "my_pca" function
	<headers> is a pd.Series containing the names of the samples in following format: "POPNAME_INDEX"
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	aa = read_my_labels(header)
	left = 0
	for i in aa.values:
		plt.plot(genotypes_PCA[left:left+i[0],0], genotypes_PCA[left:left+i[0],1], '.', label=i)
		left = left + i[0]
	plt.savefig(plot_filename)

def find_majority( mia_lista, error, exclude ):
	'''
	Finds the element with the highest frequency in <mia_lista> that is not present in <exclude>. 
	<mia_lista> is kmeans.labels_ list
	<error> is a python list that for each different label in <mia_lista> will return the percentage of individuals that do not have the major label
	<exclude> is a python list that contains already excluded labels
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	count_dict=cl.Counter( mia_lista )
	b=(max( [[count_dict[kleidi], kleidi] for kleidi in count_dict.keys() if kleidi not in exclude]))
	error.append( mia_lista.shape[0] - mia_lista[ mia_lista == b[1] ].shape[0] )
	exclude.append(b[1])
	return error, exclude

def k_means(np_table):
	'''
	Performs kmeans analysis with clusters = the number of different populations and returns the percentage of individuals labeled incorrectly (<total_error>)
	<np_table> is the output of "my_pca" function 
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	header=pd.Series(np_table[:,0])
	aa = read_my_labels(header)
	genotypes_PCA=np_table[:,1:].astype(float)
	kmeans = KMeans(n_clusters= len(aa.keys()), random_state=0).fit(genotypes_PCA)
	left = 0
	error = []
	exclude = []
	for i in aa.keys():
		error, exclude = find_majority( kmeans.labels_[aa[i]], error, exclude )
	total_error = sum(error)/kmeans.labels_.shape[0]
	return total_error
 
def find_ratio(pop_vecs, sample_table, population_list, population_sizes, independent = None, min_af = 0, iter_time=1):
	"""
	Find the minimum # of dependent SNPs needed in the sample, so that the error is less than 0.1. It outputs the # of SNPs and the eror_rate.
	If the populations do not pass the error criterion after 50 cycles, the function returns the current error rate as well.
	<pop_vectors> is the pandas DataFrame containing genotypes frequencies
	<sample_table> is a pandas DataFrame containing information from the sample information file
	<population_list> is a python list with all the names of the populations provided by the command line
	<populaiton_sizes> is a python dictionary with keys = pop_names and values = how many individuals are to be simulated for each population
	<independent> is the number of SNPs to be simulated. It is provided by command line
	<min_af> is the lnowest acceptable frequency of a SNP
	authors: Maria Malliarou and Stefanos Papadantonakis
	"""
	print ("Finding error ratio in these populations",population_list,". Iteration:",iter_time,"\n")
	
	simulation_table=simulation_independent(pop_vecs,
						sample_table,
						population_list,
						population_sizes,
						args.independent,
						min_af=args.MINIMUM_AF)
	simulation_table = np.vstack((np.array(simulation_table.columns),simulation_table.as_matrix()))
	pca_table = my_pca(simulation_table[1:,])
	k_means_value = k_means(np.column_stack((simulation_table[0,:],pca_table)))
	counter = 0
	while True:
		counter += 1
		dependent_simulation_table=simulation_dependent(pop_vecs,
								population_list,
								population_sizes,
								int( independent/10 ),
								min_af=args.MINIMUM_AF).as_matrix()
		simulation_table = np.vstack((simulation_table,dependent_simulation_table))
		pca_table = my_pca(simulation_table[1:,])
		error_rate = k_means(np.column_stack((simulation_table[0,:],pca_table)))
		print(error_rate)
		if error_rate<0.1:
			return [int(simulation_table.shape[0]-independent-1), float(error_rate) ]
		elif (counter >= 50):
			return [float(simulation_table.shape[0]-independent-1), float(error_rate) ]
		else:
			k_means_value = error_rate

def dendrogram (pop_vectors, sample_table, independent=None, min_af = 0, iterations=1):
	'''
	Constructs a dendrogram between all populations described in the <sample_table> using the error metrics produced by the "find_ratio" funtion
	<pop_vectors> is the pandas DataFrame containing genotypes frequencies
	<sample_table> is a pandas DataFrame containing information from the sample information file
	<independent> is the number of SNPs to be simulated. It is provided by command line
	<min_af> is the lnowest acceptable frequency of a SNP
	authors: Maria Malliarou and Stefanos Papadantonakis
	'''
	population_sizes=dict(sample_table["pop"].value_counts())
	populations=list(population_sizes.keys())
	#populations= ["ITU","FIN","YRI","CLM","KHV"]
	list_of_pairs = [[populations[p1], populations[p2]] for p1 in range(len(populations)) for p2 in range(p1+1,len(populations))]
	d = [ np.array([find_ratio(pop_vectors, 
						sample_table,
						pairs, 
						{pairs[0]:population_sizes[pairs[0]], pairs[1]:population_sizes[pairs[1]]},
						independent,
						min_af)
						for i in range(0,iterations)]).mean(axis = 0)
			for pairs in list_of_pairs]
	ytdist = np.array( [i[0] * i[1] for i in d ] )
	Z = hierarchy.linkage(ytdist, 'single')
	plt.figure(figsize=(15, 5))
	dn = hierarchy.dendrogram(Z, labels=populations)
	plt.show()
	return ytdist


################------- THIS IS HOW WE DO IT -----------################

action=args.action

if action == "VCF_INFO" :	###part1  
	try : vcf_info(args.vcf, start = args.START, end = args.END)
	except Exception as e:
		print('--------------------------------------------')
		print ('--vcf must be a .vcf file or a zipped .vcf.gz file')
		print('--------------------------------------------')
		raise e
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths in order to run this action.Please check for typing errors or your current folder" )
		print('--------------------------------------------')
		raise f
elif action == "VALIDATE_SAMPLE_INFO" : ####part 2
	try : validate_sample(args.vcf, args.sample_filename)
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths in order to run this action.Please check for typing errors or your current folder" )
		print('--------------------------------------------')
		raise f
elif action== "SAMPLE_INFO" :	####part 2
	try : sample_info(args.sample_filename)
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths in order to run this action.Please check for typing errors or your current folder" )
		print('--------------------------------------------')
		raise f
elif action == "SIMULATE" :  ####parts 3-5
	try:
		freq_table=read_or_create_frequency_file(args.sample_filename, args.vcf)
		populations=[i[0] for i in args.population]
		population_sizes={i[0]:int(i[1]) for i in args.population}
		pop_vectors={i:freq_table[i] for i in populations}
		sample_table=pd.read_csv(args.sample_filename, sep="\t", header = 0)
		print (pop_vectors.keys())
		if args.independent == None :
			simulation_table=simulation_dependent(pop_vectors, populations, population_sizes, args.SNPs, min_af=args.MINIMUM_AF)
			simulation_table.to_csv(args.output, sep="\t", mode='w', index=False)
		elif args.independent>0 :
			simulation_table_dependent=simulation_dependent(pop_vectors, populations, population_sizes, args.SNPs, min_af=args.MINIMUM_AF)
			simulation_table_independent=simulation_independent(pop_vectors, sample_table, populations, population_sizes, args.independent, min_af=args.MINIMUM_AF)
			pd.concat((simulation_table_dependent,simulation_table_independent), axis=0).to_csv(args.output, sep="\t", mode='w', index=False)
		else:
			raise Exception ("The argument --independent must me a positive number") 
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths  for the arguments --vcf and --sample_filename in order to run this action.Please check for typing errors or your current folder","\n" )
		print('--------------------------------------------')
		raise f
	except TypeError as t :
		print('--------------------------------------------')
		print ("Possible mistakes:,","\n","the argument --SNPs must me a positive integer","\n","the argument --population must have a specific population name and a positive integer","\n")
		print('--------------------------------------------')
		raise t
	except Exception as e :
		print('--------------------------------------------')
		print ("The population name you typed is either spelled wrong or doesnt exist in the sample_filename or your vcf file. You can run the validate info action to make sure it exists in your files")
		print('--------------------------------------------')
		raise e
elif action == "PCA" : ### part 6
	try :
		np_table=np.loadtxt(args.input_filename, delimiter="\t", dtype="object")
		header=pd.Series(np_table[0,])
		pca_table=my_pca(np_table[1:,])
		pca_plot(pca_table, header)  ### for the plot 
		header=np.array(header)
		np.savetxt(args.PCA_filename, np.column_stack((header,pca_table)) ,fmt="%s\t%f\t%f")
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths  for the argument --input_filename in order to run this action.Please check for typing errors or your current folder","\n" )
		print('--------------------------------------------')
		raise f
elif action == "CLUSTER" : ##part 7
	try:
		pca_table=np.loadtxt(args.PCA_filename, delimiter = '\t', dtype = 'object')
		print("Error rate is {}".format(k_means(pca_table)))
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths  for the argument --PCA_filename in order to run this action.Please check for typing errors or your current folder","\n" )
		print('--------------------------------------------')
		raise f
elif action == "FIND_RATIO" : ##part 8
	try:
		freq_table=read_or_create_frequency_file(args.sample_filename, args.vcf)
		populations=[i[0] for i in args.population]
		population_sizes={i[0]:int(i[1]) for i in args.population}
		pop_vectors={i:freq_table[i] for i in populations}
		sample_table=pd.read_csv(args.sample_filename, sep="\t", header = 0)
		ratio_table = np.array( [find_ratio(pop_vectors, sample_table, populations, population_sizes, args.independent, args.MINIMUM_AF, iter_time=i) 
								for i in range(0,args.iterations) ] )
		print ("SNPs needed", np.mean([i[0] for i in ratio_table]))
		if len(ratio_table) > 1:
			print( 'Confidence Interval: ')
			print (st.t.interval(0.95, float(len(ratio_table)-1), loc=np.mean([i[0] for i in ratio_table]), scale=st.sem([i[0] for i in ratio_table])))
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths  for the arguments --vcf and --sample_filename in order to run this action.Please check for typing errors or your current folder","\n" )
		print('--------------------------------------------')
		raise f
	except TypeError as t :
		print('--------------------------------------------')
		print ("Possible mistakes:,","\n","the argument --SNPs must me a positive integer","\n","the argument --population must have a specific population name and a positive integer","\n")
		print('--------------------------------------------')
		raise t
	except Exception as e :
		print('--------------------------------------------')
		print ("The population name you typed is either spelled wrong or doesnt exist in the sample_filename or your vcf file. You can run the validate info action to make sure it exists in your files")
		print('--------------------------------------------')
		raise e
elif action == "DENDROGRAM" : ##part 10
	try:	
		freq_table=read_or_create_frequency_file(args.sample_filename, args.vcf)
		sample_table=pd.read_csv(args.sample_filename, sep="\t", header = 0)
		pop_vectors=freq_table.to_dict("series")
		q=dendrogram(pop_vectors, sample_table, args.independent, args.MINIMUM_AF, iterations=args.iterations)
		np.savetxt("dendrogram_results.txt",  q ,fmt="%f")
	except FileNotFoundError as f :
		print('--------------------------------------------')
		print ("You need to put the right files with the corect paths  for the arguments --vcf and --sample_filename in order to run this action.Please check for typing errors or your current folder","\n" )
		print('--------------------------------------------')
		raise f
	except TypeError as t :
		print('--------------------------------------------')
		print ("Pssible mistakes:,","\n","the argument --SNPs must me a positive integer","\n","the argument --population must have a specific population name and a positive integer","\n")
		print('--------------------------------------------')
		raise t






