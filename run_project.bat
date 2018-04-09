::part 1
py python_project1.py --vcf 1kgp_chr1.vcf --action VCF_INFO
::part 2
py python_project1.py --vcf chr22_100000_lines.vcf --sample_filename sample_information.csv --action VALIDATE_SAMPLE_INFO
py python_project.py --sample_filename sample_information.csv --action SAMPLE_INFO
::part 3
py project.py --vcf 1kgp_chr1.vcf --sample_filename sample_information.csv --population GBR 100 --SNPs 100 --action SIMULATE --output outputfile.vcf
::part 4
py python_project1.py --vcf chr22_100000_lines.vcf --sample_filename sample_information.csv --population YRI 100 --population FIN 100 --SNPs 100 --action SIMULATE --output outputfile_test.vcf
::part 5 

::py python_project1.py --vcf chr22_100000_lines.vcf --sample_filename sample_information.csv --population GBR 100 --SNPs 100 --output outputfile.vcf -- population FIN 100 -- population YRI 200 --independent 3000 --MINIMUM_AF 0.3
py python_project1.py --vcf 1kgp_chr1.vcf --sample_filename sample_information.csv --population GWD 100 --population YRI 200 --SNPs 500 --output outputfile.vcf --action SIMULATE --independent 300 --population IBS 300 --MINIMUM_AF 0.3
::part 6
py python_project.py --input_filename outputfile.vcf --PCA_filename pca.txt --PCA_plot pca.pdf --action PCA
::part 7
py python_project.py --PCA_filename pca.txt  --action CLUSTER
::part 8
py python_project1.py --vcf chr22_100000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --action FIND_RATIO --independent 100 --MINIMUM_AF 0.1 --iterations 2
py python_project.py --vcf 1kgp_chr1.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --action FIND_RATIO --independent 100 --MINIMUM_AF 0.1 --START 16050678 --END 16206607 --iterations 10 --jaccard Tre

py python_project.py --vcf 1kgp_chr1.vcf --action "VCF_INFO"
::part 10
py python_project1.py --vcf chr22_100000_lines.vcf --sample_filename sample_information.csv --action DENDROGRAM --independent 200 --MINIMUM_AF 0.2