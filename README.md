# PGx_project

## You will need the following files:
### Files from OpenArray
1. CYP_alleles.txt: TXT file results exported from AlleleTyper
2. Genotype_Matrix.csv: Genotype matrix exported from OpenArray software

### Additional tab-separted files to be prepared
Each files contains a header and the first column containing sample ids (MUST BE the same as reported in OpenArray files)
<br>3. SLC6A4_LS.csv: patients IDs (col1) and L/S genotype (col2) for SLC6A4 gene
<br>4. 1st_drug.csv: patients IDs (col1) and first drugs used (col2)

### Resource files present in the resources folder
<br>5. AssayID_to_rsID.txt: assayID (col1), gene_name (col2) and rsID (col3)
<br>6. Metabolizer_table.csv: cyp gene (col1), star allele (col2) and metabolyzer pheno (col3)

## Steps to generate PGx report
1. Use Prepare_data.pl to prepare input table for report generation
<br>`perl Process_CYP_alleles.pl CYP_alleles.txt Genotype_Matrix.csv SLC6A4_LS.csv 1st_drug.csv resources/AssayID_to_rsID.csv resources/Metabolizer_table.csv output_file.tsv`
