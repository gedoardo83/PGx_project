# PGx_project

The provided script process input data from AlleleTyper, TaqMan Genotyper and SLC6A4 L/S variants and save a summary table of data. Web app input tableare udpated automatically. If previous data is present, the table will be backed up.

## Input files:
See example of the input files in example_input_files folder 

1. Alleletyper: TXT file results exported from AlleleTyper
2. Genotype Matrix: Genotype matrix exported from TaqMan Genotyper
3. SLC6A4: tab-separated containing patients IDs (col1) and L/S genotype (col2)
4. Sample info: tab-separated containing: clinicianID, patientID, experimental group and first drug. Columns are supposed to be in this order, but you can provide specific column indeces eventually

## Resource files in the resources folder
These files are already present in the folder and are not supposed to be modified
* AIF file: AIF file of the custom openarray design
* CYP2D6_Metabolyzer.csv: Metabolyzer phenotypes for the different CYP2D6 diplotypes
* CYP2C19_Metabolyzer.csv: Metabolyzer phenotypes for the different CYP2D6 diplotypes
* CYP2D6_AF.tsv: Allele frequencies for the CYP2D6 diplotypes
* Two Guideline files containing information for translating genotypes into reccomenations

## Usage
`python Prepare_data.py --sample_data Sample_info.csv --genos TaqMan_genotyper_results.txt --cyp_alleles Allele_typer_result.txt --slc6a4 SLC6A4_LS_alleles.txt`
 
### Mandatory arguments
* --sample_data: Sample data including sampleID,clinicianID,first drug used,exp group
* --genos: Genotype matrix for SNP from taqman genotyper
* --cyp_alleles: CYP alleles from Allele Typer
* --slc6a4: Long / short alleles for SLC6A4

### Optional arguments
* --sample_cols: Comma-separated col numbers for ClinicianID, SampleID, PGx_group and 1stdrug columns, in this order
* --resources: Resources folder for metabolyzer table, AIF file and haplotype frequencies
* --ntc_id: Sample ID for negative control, this will be removed from tables

### Customization
All folders locations and resource files names can be modifeid as needed editing the configuration section at the beginning of the Prepare_data.py script

## Output
* data/PGx_data.tsv updated: contains the processed dataset for the current input data (sample information, cyp alleles, cyp metabolyzer phenotypes, additional genotypes)
* web_app tables updated 
* your input files moved to input_data folder with the current date prefix for archiving

## Restore the last version of data
During process, the last web app datasets are saved in app_backup folder with the prefix of the current time.
If anything go wrong processing new files, you MUST replace the newly created .csv files with the old ones. This will ensure that the web app still run while you fix the problem.
