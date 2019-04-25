# PGx_project

Do not alter directory structure and do not change filenames for files in resources and web_app folder

## You will need the following files. Filenames are not important:
### Files from OpenArray
1. Alleletyper: TXT file results exported from AlleleTyper
2. GenoMatrix: Genotype matrix exported from OpenArray software

### Additional tab-separted files to be prepared
Each files contains a header and the first column containing sample ids (MUST BE the same as reported in OpenArray files)
3. SLC6A4: patients IDs (col1) and L/S genotype (col2) for SLC6A4 gene
4. 1stDrug: patients IDs (col1) and first drugs used (col2)

### Resource files present in the resources folder
AssayID_to_rsID.txt: assayID (col1), gene_name (col2) and rsID (col3)
<br>Metabolizer_table.csv: cyp gene (col1), star allele (col2) and metabolyzer pheno (col3)
<br>Two Guideline files containing information for translating genotypes into reccomenations

## Steps to update PGx report web app
1. Update the samples_info.csv file
This is a 3-columns tab-separated file containing subject informations and header. Columns names must be as follows: Clinician, Sample, Group_PGx.
<br>Clinician = Clinician unique code
<br>Sample = Sample unique code (MUST BE the same used in genotyping)
<br>Group_PGx = 1 for having a report, 0 for exluded from report

2. Use the Process_results.sh script
The script take your genotyping results as input and create tables need for the PGx report. New results are added to previous dataset eventually present in the folders and everything is automatically updated and deployed to web application
<br>`bash Process_results.sh  AlleleTyper_file GenoMatrix_file SLC6A4_file 1stDrug_file`

At the end the overall dataset of relevant genotypes is created in sample_tables/ALL_samples_genos.csv
