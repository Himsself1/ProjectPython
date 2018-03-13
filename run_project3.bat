::part 1
py project3.py --vcf chr22_3000_lines.vcf --action VCF_INFO
::part 2
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --action VALIDATE_SAMPLE_INFO
py project3.py --sample_filename sample_information.csv --action SAMPLE_INFO
::part 3
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population GBR 100 --SNPs 100 --action SIMULATE --output outputfile.vcf
::part 4
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population GBR 100 --population FIN 100 --SNPs 100 --action SIMULATE --output outputfile.vcf
::part 5 

::py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population GBR 100 --SNPs 100 --output outputfile.vcf -- population FIN 100 -- population YRI 200 --independent 3000
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --SNPs 500 --output outputfile.vcf --action SIMULATE --independent 300 --population GBR 300
::part 6
py project3.py --input_filename outputfile.vcf --PCA_filename pca.txt --PCA_plot pca.pdf --action PCA
::part 7
py project3.py --PCA_filename pca.txt  --action CLUSTER
::part 8
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --action FIND_RATIO --independent 100 --MINIMUM_AF 0.1 --START 16050678 --END 16206607 --iterations 10
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --action FIND_RATIO --independent 100 --MINIMUM_AF 0.1 --START 16050678 --END 16206607 --iterations 10 --jaccard Tre

py project3.py --vcf chr22_3000_lines.vcf --action "VCF_INFO"
::part 10
py project3.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --action DENDROGRAM 