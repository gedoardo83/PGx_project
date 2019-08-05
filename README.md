# PGx_project

Do not alter directory structure and do not change filenames for files in resources and web_app folder

## You will need the following files. Filenames are not important:
### Files from OpenArray
1. Alleletyper: TXT file results exported from AlleleTyper.

This file needs to be processed to remove ^M character. Use the following command

`sed -i "s/^M//g" Alleletyper_file`

To enter ^M, type CTRL-V, then CTRL-M. That is, hold down the CTRL key then press V and M in succession.
2. GenoMatrix: Genotype matrix exported from TaqMan Genotyper

### Additional tab-separted files to be prepared
Each files contains a header and the first column containing sample ids (MUST BE the same as reported in OpenArray files)

3. SLC6A4: patients IDs (col1) and L/S genotype (col2) for SLC6A4 gene
4. 1stDrug: patients IDs (col1) and first drugs used (col2)

<br>
Examples of all input files are provided in example_input_files folder

## Resource files present in the resources folder
These files are already present in the folder and are not supposed to be modified
* CustomArray_AIF.txt: AIF file of the custom openarray design
* Metabolizer_table.csv: cyp gene (col1), star allele (col2) and metabolyzer pheno (col3)
* Two Guideline files containing information for translating genotypes into reccomenations

## Procedure to generate new reports
### 1. Update the samples_info.csv file present in the main folder
This is a 3-columns tab-separated file containing subject informations and header. Columns names must be as follows: Clinician, Sample, Group_PGx.
* Clinician = Clinician unique code
* Sample = Sample unique code (MUST BE the same used in genotyping)
* Group_PGx = 1 for having a report, 0 for excluded from report

### 2. Use the Process_results.sh script
The script takes your genotyping results as input and create tables needed for the PGx report. New results are added to previous datasets eventually present in the folders and everything is automatically updated and deployed to web application
`bash Process_results.sh AlleleTyper_file TaqManGenotyper_file SLC6A4_file 1stDrug_file`

<br>

At the end pf process you will have
* sample_tables/[presentdate]_samples_genos.csv: contains the processed dataset for the current input data
* sample_tables/ALL_samples_genos.csv: contains the overall dataset used for the report for all samples analyzed so far (genotypes, CYP alleles, 1st drug)
* your input files are moved to input_data folder with the current date prefix for archiving

## Restore the last version of data
During process, the last datasets are stored as .old files in 
* samples_tables folder: samples_tables/ALL_samples_genos.old
* web_app folder: web_app/samples_genos.old

If anything go wrong during process of new files, you MUST replace the newly created .csv files with the .old ones. This will ensure that the web app still run while you fix the problem. DO NOT attempt 2 subsequent updates that are likely to fail before having restored the backup files.
