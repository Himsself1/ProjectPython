###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389

!pip install numpy

#------------our imports----------------------

import argparse
import gzip
import re
import numpy as np
import assert



'lets\ttest\tthis\tshit\nthis\tis\tan\t\array'



#------------our arguments-----------------------
parser = argparse.ArgumentParser()

#part1
parser.add_argument('--vcf', nargs=1)  ###θα παίρνει ένα vcf ή ένα vcf.gz αρχείο
parser.add_argument('--action', nargs=1)  ##Το πρόγραμμά σας θα τυπώνει πόσα SNPs έχει το vcf και πόσα samples.



##example > python project.py --vcf 1kgp_chr1.vcf.gz --action VCF_INFO

#File has 2502 samples
#File has 1234567 SNPs

#part2
parser.add_argument('--sample_filename', nargs=1)

##--action == SAMPLE_INFO , ==VALIDATE_SAMPLE_INFO





##part3

parser.add_argument('--population', nargs=2)  ##POPULATION_NAME, NUMBER_OF_SAMPLES
parser.add_argument('--SNPs', nargs=1)
parser.add_argument('--output', nargs=1)
### ++++ action= SIMULATION
## > python project.py --vcf 1kgp_chr1.vcf.gz --sample_filename sample_information.csv --population GBR 100 --SNPs 1000 --output random_genotypes.vcf --action SIMULATE

##part4
###multiple population samples

##part5
parser.add_argument('--independent', nargs=1)  #<NUMBER_OF_INDEPENDENT_SNPs bases in all populations given

##part6----PCA!!

parser.add_argument('--input_filename', nargs=1)
parser.add_argument('--PCA_filename', nargs=1)
parser.add_argument('--PCA_plot', nargs=1)
##action = PCA


##part7   k-means

###--PCA_filename , --action = CLUSTER

##part8
#--independent , --action FIND_RATIO 

##part9
parser.add_argument('--iterations', nargs=1)  #NUMER_OF_ITERATIONS
parser.add_argument('--MINIMUM_AF', nargs=1)  #MINIMUM_ALLELE_FREQUENCY


parser.add_argument('----START', nargs=1)
parser.add_argument('--END', nargs=1)

#part10
#--action DENDROGRAM

#part11 (optional)  jaccard index

input_file = '1kgp_chr1.vcf.gz'

with gzip.open(input_file, 'rt') as arxeio:
    print( arxeio.readline() )


## Downloading files

args = parser.parse_args()
input_file = args.vcf

############# Reading the file #######################
if re.match('.+vcf.gz$', input_file):
    with gzip.open(input_file, 'rt') as arxeio:
        mia_metablhth = diavazw_vcf( arxeio )
elif re.match('.+.vcf$', input_file):
    with open(input_file, 'rt') as arxeio:
        mia_metablhth = diavazw_vcf( arxeio )
else:
    print( 'WTF IS THIS FILE: {}'.format(input_file) )
    assert(0)
#####################################################

print (sample_data)

########################## TEST ####################################

############# PART 1 #####################
input_file = 'chr20_lines.vcf'
file_name = open( input_file, 'rt' )
header, data = diavazw_vcf(file_name)
info_vcf( data )
file_name.close()
#########################################

############## PART 2 ###################




sample_filename = 'sample_information.csv'
sample_file = open(sample_filename)
sample_data = sample_file.read()
sample_file.close()
sample_array=np.loadtxt("sample_iformation.csv", delimiter="\t", skiprows = 1)

sample_data_splitted = [x.split() for x in sample_data.split('\n')]

print(sample_data_splitted)

setted_population=set(sample_array['super_pop'])
print(setted_population) 



print(header)

###################################################################

def diavazw_vcf( file_name ):
    '''
    Oti leei to onoma
    '''
    ## Agnooume tis grammes pou arxizoun me '##'
    data = file_name.readline()
    while data[:2] == '##':
        data = file_name.readline()
    ## Ka8arizoume to arxeio
    header = data
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




    
    
    # for i in data :
    #     info_dict = dict([x.split('=') for x in info.split(';')])
    #     if x[] :
    #         counter++

    
