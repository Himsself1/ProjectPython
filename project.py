###reference material: https://gist.github.com/kantale/81d7d728c22fb35d77112c3633e17389



#------------our imports----------------------

import argparse
import gzip




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



args = parser.parse_args()




print (sample_data)
