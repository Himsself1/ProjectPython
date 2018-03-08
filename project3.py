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
from sklearn.metrics.pairwise import pairwise_distances
from scipy.cluster import hierarchy
import scipy.stats as st
import math as math
import itertools
##################

##################### ------------our arguments-----------------------###############
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', type=str, help="Type your .vcf or .vcf.gz file. If your file is in different destination,you need to type the path too")
parser.add_argument('--action', help="Type what action you want this program to perform", choices=["VCF_INFO","SAMPLE_INFO","VALIDATE_SAMPLE_INFO","SIMULATE","PCA","CLUSTER","FIND_RATIO","DENDROGRAM"])
parser.add_argument('--sample_filename', type=str, help="Type the sample file")
parser.add_argument('--population', nargs = 2, action='append', help="Type the name of the population and the desired sample number")
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

def read_vcf_file(file_name, start=0, end=None):
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
	if end != None and end > start : 
		return df[(~df['INFO'].str.contains("MULTI")) &
				(df['POS'] > start) & 
				(df['POS'] < end) ]  ###dataframe contains only SNP's
	elif end == None :
		return df[(~df['INFO'].str.contains("MULTI")) &
		(df['POS'] > start) ]  ###dataframe contains only SNP's
	else :
		raise Exception("Invalid File Extension. File needs to end with .vcf or .vcf.gz")

def vcf_info(file_name, start = 0, end = None) :
	df=read_vcf_file(file_name, start, end)
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
	z = np.array(z)
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 2
	return z

def convert_to_jaccard (z) :
	z[ z == '0|0' ] = 0
	z[ z == '1|0' ] = 1
	z[ z == '0|1' ] = 1
	z[ z == '1|1' ] = 1
	return z

def jaccard_index (table) :
	jac_table=convert_to_jaccard(table)
	return 1 - pairwise_distances(table.T, metric = "hamming")

def find_non_zero_AFs( dict_of_pop_freqs, min_af = 0.0 ):
	'''
	Returns the indexes of acceptable SNPs
	<min_af> is the lowest acceptable frequency
	<dict_of_pop_freqs> is a dictionary with keys = population names and values = frequency vectors for SNPs
	values are numpy arrays
	'''
	table = np.array([ i for i in dict_of_pop_freqs.values() ]).sum( axis = 0 )/len( dict_of_pop_freqs.keys() )
	return np.array( np.where( table >= min_af )[0] ) 

def create_genotype( frequenc, sample_number ):
	'''
	Generates a genotype of size <sample_number> with snp frequency = <frequenc>
	<frequenc> must be a float
	<sample_number> must me float converible
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
	<data> is a pandas Data Frame
	<individuals> is a 1D pandas Data Frame with the names of individuals of the <data> 
	'''
	temp_data = convert_genotype_to_number(data[individuals])
	frequency = temp_data.sum( axis = 1 )/float(2*temp_data.shape[1])
	return frequency

def population_columns( sample_matrix, pop_name ):
	'''
	Returns the names of all individuals of population <pop_name>
	'''
	if sample_matrix['pop'].str.contains(pop_name).sum() == 0:
		raise Exception ('Error, Wrong population name: {}'.format(pop_name))
	else:
		return sample_matrix[sample_matrix['pop'].str.contains(pop_name)]['sample']


data=read_vcf_file(args.vcf, start = args.START, end = args.END)
df = pd.read_csv('sample_information.csv', sep="\t", header = 0)
populations=["GBR","FIN","YRI"]
def simulation_dependent( data, sample_matrix, pop_names, sample_sizes, total_snps, min_af = 0 ):
	'''
	Simulates dependent genotypes
	<data> is the pandas DataFrame containing genotypes
	<sample_matrix> is a pandas DataFrame containing information from the sample information file
	<pop_names> is a python list with all the names of the populations provided by the command line
	<sample_sizes> is a python dictionary with keys = pop_names and values = how many individuals are to be simulated
	for each population
	<total_snps> is an integer containing the number of snps to be simulated, provided by command line
	<min_af> is the lowest acceptable frequency of a SNP
	'''
	pop_vecs={}
	for i in pop_names:
		noumera = population_columns( df, i )
		pop_vecs[i]=calculate_frequencies( data, noumera )
	index=find_non_zero_AFs(pop_vecs, min_af)
	q=np.random.choice(index, size = total_snps, replace = True)
	output_dict = { i: pd.DataFrame([create_genotype(pop_vecs[i][ind], sample_sizes[i]) for ind in q],
									columns=[str(pop_name)+"_"+str(k) 
									for k in range(1,sample_sizes[i]+1)]) 
					for i in pop_names }
	finale=pd.concat(output_dict.values(), axis=1)
	return finale

sample_size={k : 100 for k in populations}

## Concatenating all the genotypes

