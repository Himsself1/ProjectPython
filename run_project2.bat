::py project2.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population GBR 100 --SNPs 100 --output outputfile.vcf -- population FIN 100 -- population YRI 200 --independent 3000
py project2.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --SNPs 500 --output outputfile.vcf --action SIMULATE --independent 3000
::part 6
py project2.py --input_filename outputfile.vcf --PCA_filename pca.txt --PCA_plot pca.pdf --action PCA
::part 7
py project2.py --PCA_filename pca.txt  --action CLUSTER
::part 8
py project2.py --vcf chr22_3000_lines.vcf --sample_filename sample_information.csv --population FIN 100 --population YRI 200 --action FIND_RATIO --independent 100
