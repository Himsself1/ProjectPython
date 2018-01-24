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
## ------------our arguments-----------------------
parser = argparse.ArgumentParser()

parser.add_argument('--vcf', nargs = 1)  ###θα παίρνει ένα vcf ή ένα vcf.gz αρχείο
parser.add_argument('--action', nargs = 1)  ##Το πρόγραμμά σας θα τυπώνει πόσα SNPs έχει το vcf και πόσα samples.


args = parser.parse_args()

input_file = 'chr22_25_lines.vcf'
## input_file = args.vcf


def read_vcf_file(file_name):
    '''
    Create a string with the header and a list with the data that only contains SNPs
    '''
    ## Agnooume tis grammes pou arxizoun me '##'
    data = file_name.readline()
    while data[:2] == '##':
        data = file_name.readline()
    ## Ka8arizoume to arxeio
    header=data.rstrip("\n").split("\t")
    data = file_name.read()
    data = data.split('\n')
    data_clean = [ i.split('\t') for i in data if 'VT=SNP' in i ]
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

    
print(data)

print(len(data))
    
#################### TEST #############################



#######################################################

#################### Part1 ############################
info_vcf( data ) ## 8elei akoma na to kanoume me argv
#######################################################

###################--Part2--##################################
##sample file
sample_filename=args.sample_filename()
sample_filename = 'sample_information.csv'

def sample_info(sample_file) :
    '''
    Returns information of the sample file
    '''
    sample_array = np.loadtxt(sample_file, delimiter = "\t", skiprows = 1, dtype = "object" )    
    unique_super_pop = list(set(sample_array[:,2]))
    unique_pop = list(set(sample_array[:,1]))
    population = Counter(sample_array[:,1])
    super_populations = Counter(sample_array[:,2])
    which_pop = dict(zip(sample_array[:,1], sample_array[:,2]))
    print ("File has", len(super_populations.keys()), "Areas.")
    for i in range(0,len(unique_super_pop)) :
        print ("Area", i+1, "is",
               unique_super_pop[i],
               "and contains",
               super_populations[unique_super_pop[i]],
               "samples splitted in the following populations:")
        for d in unique_pop:
            if which_pop[d] == unique_super_pop[i]:
                print (d,population[d], "samples")
                
def validate_sample(sample_file, sample_ids):
    '''
    Check if sample file has the same samples with vcf file
    '''
    sample_array = np.loadtxt(sample_file, delimiter="\t", skiprows = 1, dtype="object" )
    samples = sample_array[:,0]
    id_array = np.array(sample_ids)
    if np.array_equiv(np.sort(samples),np.sort(id_array)):
        print ("Everything is OK!")
    else :
        b = np.setdiff1d(id_array,samples)
        q = np.setdiff1d(samples,id_array)
        if len(b) != 0:
            print( "These samples are present in VCF but not in SAMPLE:", b )
        if len(q) != 0:
            print( "These samples are present in SAMPLE but not in VCF:", q )



###################--Part2--##################################

## sample_filename=args.sample_filename()


################## Part3 ##########################


parser.add_argument('--population', nargs = 2)  ##Το πρόγραμμά σας θα τυπώνει πόσα SNPs έχει το vcf και πόσα samples.
parser.add_argument('--SNPs', nargs = 1) 
parser.add_argument('--output', nargs = 1)


def generate_pop_freq( data, header, pop_info, pop_name ):
    '''
    Generates a frequency vector for each SNP for individuals in a given population
    '''
    my_data = pd.DataFrame(data, columns = header)
    ## Stores the decription of individuals
    pop_parameters = np.loadtxt( pop_info, delimiter = "\t", skiprows = 1, dtype = "object" ) 
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


def writeSNPs( my_sz, my_snps, frequency, my_output_file ):
    '''
    Generates the file that stores the randomized genotype samples
    '''
    output = []
    for it in range(0,my_snps):
        temp = random.choice(frequency)
        ones = int(temp*2*my_sz // 1) ## Posa 1 8a uparxoun sto dataset? Takes the integer part of the number
        malakia1 = [0]*2*my_sz
        a = random.sample( range(2*my_sz), ones )
        for i in a:
            malakia[i] = 1 ## Ftiaxnw tous assous
        genotype = []
        for k in range(0,len(malakia1)-1,2):
            genotype.append(str(malakia[k])+'|'+str(malakia[k+1])) ## Kanw concatenate to genotype
        output.append(genotype)
    with open(my_output_file, 'w') as f :
        ## NO HEADER AND CAPES!!!!
        for i in output:
            f.write("\t".join(i)+"\n")


#### TEST ######

my_population = "MSL"
my_no_samples = 100
my_snps = 1000
my_output = "random_genotypes.vcf"

sample_filename = 'sample_information.csv'
sample_info( sample_filename )

test = generate_pop_freq( data, header, sample_filename, my_population )
print(test)

################


############ Part 4 ###############
