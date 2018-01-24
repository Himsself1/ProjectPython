#Copyright (C) 2017-2018 Alexandros Kanterakis, mail:kantale@ics.forth.gr
###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389


#########################----IMPORTS-------------------------------

import argparse
import gzip
import re
import numpy as np
from collections import Counter
#------------our arguments-----------------------
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', nargs=1)  ###θα παίρνει ένα vcf ή ένα vcf.gz αρχείο
parser.add_argument('--action', nargs=1)  ##Το πρόγραμμά σας θα τυπώνει πόσα SNPs έχει το vcf και πόσα samples.


args = parser.parse_args()


input_file = 'chr22_25_lines.vcf'
## input_file = args.vcf


def read_vcf_file(file_name):
        '''
        Create a string with the header and a list with the data
        '''
        ## Agnooume tis grammes pou arxizoun me '##'
        data = file_name.readline()
        while data[:2] == '##':
                data = file_name.readline()
	## Ka8arizoume to arxeio
        header = data.rstrip("\n").split("\t")
        print(header)
        data = file_name.read()        
	data = data.split('\n')
        data_clean = [ i.split('\t') for i in data ]
        if len(data_clean[-1]) == 0:
                data_clean = data_clean[0:len(data_clean)-1]
        return header, data_clean

def info_vcf( data ):
	'''
	Briskw posa samples kai posa SNPs exei to arxeio 
	'''
	###information = [ i[7].split(';') for i in data ]
	AC = []
	for i in data:
		info_dict = dict([j.split('=') for j in i[7].split(';')])
		if info_dict["VT"]=="SNP" :
			AC.append(int(info_dict['AC']))
	# print (info_dict)
	print( 'Minor alleles are: {}'.format(sum(AC)) )
	print( 'Samples are: {}'.format(info_dict['NS']) )


###open the file

if re.match('.+vcf.gz$', input_file):
	with gzip.open(input_file, 'rt') as arxeio:
		header, data = read_vcf_file( arxeio )
elif re.match('.+.vcf$', input_file):
	with open(input_file, 'rt') as arxeio:
		header, data = read_vcf_file( arxeio )
else :
	print ("Incorrect file format")

sample_ids = header[9:]


###################--Part2--##################################
##sample file 
#sample_filename=args.sample_filename()
sample_filename = 'sample_information.csv'

def sample_info(sample_file) :
	'''
	Returns information of the sample file
	'''
	sample_array=np.loadtxt(sample_file, delimiter="\t", skiprows = 1, dtype="object" )

	unique_super_pop=list(set(sample_array[:,2]))
	unique_pop=list(set(sample_array[:,1]))
	population=Counter(sample_array[:,1])
	super_populations=Counter(sample_array[:,2])
	which_pop=dict(zip(sample_array[:,1], sample_array[:,2]))
	print ("File has", len(super_populations.keys()),"Areas.")
	for i in range(0,len(unique_super_pop)) :
		print ("Area",i+1,"is", unique_super_pop[i],"and contains", super_populations[unique_super_pop[i]], "samples splitted in the following populations:")
		for d in unique_pop:
			if which_pop[d]==unique_super_pop[i] :
				print (d,population[d],"samples")

def validate_sample(sample_file):
	'''
	Check if sample file has the same samples with vcf file
	'''
	sample_array=np.loadtxt(sample_file, delimiter="\t", skiprows = 1, dtype="object" )
	samples=sample_array[:,0]
	id_array=np.array(sample_ids)
	if np.array_equiv(np.sort(samples),np.sort(id_array)) :
		print ("Everything is OK!")
	else :
		b=np.setdiff1d(id_array,samples)
		q=np.setdiff1d(samples,id_array)
		if len(b)!=0:
			print ("These samples are present in VCF but not in SAMPLE:",b)
		if len(q)!=0:
			("These samples are present in SAMPLE but not in VCF:",q)

			
###part1


###################--Part2--##################################
validate_sample(sample_filename)







