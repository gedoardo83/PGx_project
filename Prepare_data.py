'''
PGx report app
Author: Edoardo Giacopuzzi

This script processes input data and generate tables for the reporting app
Input are:
- sample definitions
- allele typer results for CYP2D6 and CYP2C19
- genotype matrix for SNPs
- genotype results for SLC6A4
- metabolyzer pheno table
- AIF file for the SNP panel
'''

from datetime import datetime
import argparse
import pandas as pd
import gzip, re, sys, os, itertools
from shutil import copyfile

#######################
###  Configuration  ###
#######################

RESOURCES_DIR="resources"
APP_DIR="web_app"
DATA_DIR="data"
INPUT_DIR="input_data"
BACKUP_DIR="app_backup"
LOG_DIR="logs"

CYP2D6_PHENO="CYP2D6_Metabolyzer.csv"
CYP2C19_PHENO="CYP2C19_Metabolyzer.csv"
AIF="CustomArray_AIF.txt"
CYP2D6_AF="CYP2D6_AF.tsv"

POP = 'European'
SNPS = ["rs489693","rs4713916","rs7997012","rs6295"]

now = datetime.now()
current_time = now.strftime("%Y%m%d_%H%M%S")

def tokenize(line,sep="\t"):
	line = line.rstrip('\n')
	line = line.split(sep)
	return line

#simple reader for a file plan text or gzip
def reader(file, sep, header=True, skipper=None):
	if file.endswith('gz'):
		f = gzip.open(file, 'rt')
	else:
		f = open(file)
	
	line = f.readline()
	if isinstance(skipper, int):
		for _ in range(skipper):
			line = f.readline()
	else:
		while skipper:
			line = f.readline()
	
	if header==True:
		colnames = tokenize(line, sep)
		line = f.readline()

	while line:
		line = tokenize(line, sep)
		if header == True:
			line = dict(zip(colnames, line))
		yield line
		line = f.readline()
	f.close()

#Command line arguments
parser = argparse.ArgumentParser(description='Process data for PGx report and create table used by the report app')
parser.add_argument("-s", "--sample_data", help="Sample data including sampleID,clinicianID,first drug used,exp group", action="store", required=True)
parser.add_argument("-p", "--sample_cols", help="Comma-separated col numbers for ClinicianID, SampleID, PGx_group and 1stdrug columns, in this order", action="store", default='0,1,2,3', required=False)
parser.add_argument("-g", "--genos", help="Genotype matrix for SNP from taqman genotyper", action="store", required=True)
parser.add_argument("-c", "--cyp_alleles", help="CYP alleles from Allele Typer", action="store", required=True)
parser.add_argument("-l", "--slc6a4", help="Long / short alleles for SLC6A4", action="store", required=True)
parser.add_argument("-r", "--resources", help="Resources folder for metabolyzer table, AIF file and haplotype frequencies", action="store", default=RESOURCES_DIR, required=False)
parser.add_argument("-r", "--ntc_id", help="Sample ID for negative control, this will be removed from tables", action="store", default='NTC', required=False)
args = parser.parse_args()

cyp2d6_pheno_file = RESOURCES_DIR + "/" + CYP2D6_PHENO
cyp2c19_pheno_file = RESOURCES_DIR + "/" + CYP2C19_PHENO
cyp2d6_af_file = args.resources + "/" + CYP2D6_AF
aif_file = args.resources + "/" + AIF
ntc_id = args.ntc_id

sample_cols_names = ['Clinician','Sample','Group_PGx','1st_drug']
sample_cols = args.sample_cols.split(',')

###############################
###  Read resources tables  ###
###############################

#Store metabolyzer data
print("Reading metabolyzer tables")
print("CYP2D6:", cyp2d6_pheno_file)
print("CYP2C19:", cyp2c19_pheno_file)
CYP2D6_pheno = pd.read_csv(cyp2d6_pheno_file, sep="\t", header=0)
CYP2C19_pheno = pd.read_csv(cyp2c19_pheno_file, sep="\t", header=0)

#Store AIF file
print("Reading AIF file:", aif_file)
AIF_table = pd.read_csv(aif_file, sep="\t", header=0, comment='*', encoding = "ISO-8859-1", usecols=['Assay ID','NCBI SNP Reference'])

#Store sample infos
sample_cols = [c-1 for c in args.sample_cols.split(',')]
samples_info = pd.read_csv(args.sample_data, sep="\t", header=0)
samples_info = samples_info.iloc[:,args.sample_cols.split(',')]
samples_info.columns=['Clinician','Sample','Group_PGx','first_drug']
samples_info.set_index('Sample',inplace=True)

#Read CYP2D6 AF
cyp2d6_af_table = pd.read_csv(cyp2d6_af_file, sep="\t", header=0)
cyp2d6_af_table.fillna(0, inplace=True)
cyp2d6_af_table.sort_values(POP, ascending=False, inplace=True)

#######################################################
###  Read actual results for CYPs, SNPs and SLC6A4  ###
#######################################################

#Read CYPs alleles from allele typer
CYP_alleles_table = pd.read_csv("detail_result.txt", sep="\t", skiprows=10, encoding = "ISO-8859-1", usecols=['sample ID','CYP2D6','CYP2C19'])
CYP_alleles_table.rename(columns={'sample ID' : 'Sample'}, inplace=True)
CYP_alleles_table.set_index('Sample', inplace=True)
if ntc_id in CYP_alleles_table.index.values:
	CYP_alleles_table.drop(ntc_id, inplace=True)
if 'Run date' in CYP_alleles_table.index.values:
    CYP_alleles_table.drop('Run date', inplace=True)

#Extract single diplotypes and select the one with the high mean AF
col_idx = CYP_alleles_table.columns.get_loc('CYP2D6')
for idx, value in enumerate(CYP_alleles_table['CYP2D6'].values):
    genos = value.split(', ')
    af_mean = []
    for g in genos:
        g = re.sub(r'x[0-9]+','xN', g)
        alleles = g.split('/')
        af_sum = sum(a for a in cyp2d6_af_table.loc[alleles,POP])
        af_mean.append(af_sum / 2)
    
    sort_idx = sorted(range(len(af_mean)), key=lambda k: af_mean[k])
    print(CYP_alleles_table.iloc[idx, col_idx], genos[sort_idx[0]])
    CYP_alleles_table.iloc[idx,col_idx] = genos[sort_idx[0]]

#Translate diplotypes to metabolizer phenos
CYP_alleles_table['CYP2D6_pheno'] = ""
CYP_alleles_table['CYP2C19_pheno'] = ""

cnv = re.compile('(\*[0-9]+)x([0-9]+)')
col_idx = CYP_alleles_table.columns.get_loc('CYP2D6_pheno')

missed_diplotypes = []

#CYP2D6 pheno
for idx, geno in enumerate(CYP_alleles_table['CYP2D6'].values):
    m = cnv.search(geno)
    if m is not None:
        if m.group(2) == "2":
            if geno in CYP2D6_pheno['Diplotype'].values:
                CYP_alleles_table.iloc[idx,col_idx]=CYP2D6_pheno.loc[CYP2D6_pheno['Diplotype']==geno, 'Phenotype'].values[0]
            else:
                geno_new = re.sub(r'x[3-9]+','xN', geno)
                if geno_new in CYP2D6_pheno['Diplotype'].values:
                    CYP_alleles_table.iloc[idx,col_idx]=CYP2D6_pheno.loc[CYP2D6_pheno['Diplotype']==geno_new, 'Phenotype'].values[0]
                else:
                    print("CYP2D6 dyplotype not found:", geno)
                    missed_diplotypes.append(geno)
                    CYP_alleles_table.iloc[idx,col_idx]="UNDETERMINED"
        elif m.group(2) > 2:
            geno_new = re.sub(r'x[3-9]+','xM', geno)
            if geno_new in CYP2D6_pheno['Diplotype'].values:
                CYP_alleles_table.iloc[idx,col_idx]=CYP2D6_pheno.loc[CYP2D6_pheno['Diplotype']==geno_new, 'Phenotype'].values[0]
            else:
                geno_new = re.sub(r'x[3-9]+','xN', geno)
                if geno_new in CYP2D6_pheno['Diplotype'].values:
                    CYP_alleles_table.iloc[idx,col_idx]=CYP2D6_pheno.loc[CYP2D6_pheno['Diplotype']==geno_new, 'Phenotype'].values[0]
                else:
                    print("CYP2D6 dyplotype not found:", geno)
                    missed_diplotypes.append(geno)
                    CYP_alleles_table.iloc[idx,col_idx]="UNDETERMINED"
    else:
        if geno in CYP2D6_pheno['Diplotype'].values:
            CYP_alleles_table.iloc[idx,col_idx]=CYP2D6_pheno.loc[CYP2D6_pheno['Diplotype']==geno, 'Phenotype'].values[0]
        else:
            print("CYP2D6 dyplotype not found:", geno)
            missed_diplotypes.append(geno)
            CYP_alleles_table.iloc[idx,col_idx]="UNDETERMINED"

#CYP2C19 pheno
col_idx = CYP_alleles_table.columns.get_loc('CYP2C19_pheno')
for idx, geno in enumerate(CYP_alleles_table['CYP2C19'].values):
    if geno in CYP2C19_pheno['Diplotype'].values:
        CYP_alleles_table.iloc[idx,col_idx]=CYP2C19_pheno.loc[CYP2C19_pheno['Diplotype']==geno, 'Phenotype'].values[0]
    else:
        CYP_alleles_table.iloc[idx,col_idx]="UNDETERMINED"
CYP_alleles_table.head()


#Read genotypes from TaqMan genotyper
geno_table = pd.read_csv(args.genos, sep=",", skiprows=15, encoding = "ISO-8859-1")
geno_table.rename(columns=AIF_table.to_dict()['NCBI SNP Reference'], inplace=True)
geno_table = geno_table[['Sample/Assay'] + SNPS]
geno_table.rename(columns={'Sample/Assay':'Sample'}, inplace=True)
geno_table.set_index('Sample', inplace=True)
if ntc_id in geno_table.index.values:
	geno_table.drop(ntc_id, inplace=True)
if 'Run date' in geno_table.index.values:
    geno_table.drop('Run date', inplace=True)

#Read SLC6A4 LS results
SLC6A4_table = pd.read_csv("example_input_files/SLC6A4_LS_example.tsv", sep="\t", header=0, index_col='Sample')

#Check that samples provided in results tables are in sample_info table
result = all(elem in samples_info.index.values for elem in geno_table.index.values)
if not result:
    print("Genotype table contains samples not found in sample information table")
    
result =  all(elem in samples_info.index.values for elem in CYP_alleles_table.index.values)
if not result:
    print("Allele-typer table contains samples not found in sample information table")

result =  all(elem in samples_info.index.values for elem in SLC6A4_table.index.values)
if not result:
    print("SLC6A4 table contains samples not found in sample inforamtion table")

#Merge samples_info with genotypes and CYP alleles
samples_info = samples_info.merge(CYP_alleles_table, on='Sample')
samples_info = samples_info.merge(SLC6A4_table, on='Sample')
samples_info = samples_info.merge(geno_table, on='Sample')

#Update tables for web app
app_sample_table = APP_DIR+"/samples_info.tsv"
app_geno_table = APP_DIR+"/samples_genos.tsv"
if os.path.isfile(app_sample_table):
	copyfile(app_sample_table, BACKUP_DIR+"/"+current_time+"samples_info.tsv")
	copyfile(app_geno_table, BACKUP_DIR+"/"+current_time+"samples_genos.tsv")
	samples_info.to_csv(app_sample_table,sep="\t",mode="a", header=False, columns=['Clinician','Group_PGx','first_drug'])
	samples_info.to_csv(app_geno_table,sep="\t",mode="a", header=False,columns=['CYP2D6','CYP2D6_pheno','CYP2C19','CYP2C19_pheno','SLC6A4']+SNPS)
else:
	samples_info.to_csv(app_sample_table,sep="\t", header=True, columns=['Clinician','Group_PGx','first_drug'])
	samples_info.to_csv(app_geno_table,sep="\t",header=True,columns=['CYP2D6','CYP2D6_pheno','CYP2C19','CYP2C19_pheno','SLC6A4']+SNPS)

#Update data table, store old data and save input tables
out_table = DATA_DIR+"/PGx_data.tsv"
if os.path.isfile(out_table):
	copyfile(out_table, DATA_DIR+"/"+current_time+"_PGx_data.tsv")	
	samples_info.to_csv(out_table,sep="\t",mode="a", header=False)
else:
	samples_info.to_csv(out_table,sep="\t",header=True)

copyfile(args.genos, INPUT_DIR+"/"+current_time+"_genos.txt")
copyfile(args.cyp_alleles, INPUT_DIR+"/"+current_time+"_cyp_alleles.txt")
copyfile(args.slc6a4), INPUT_DIR+"/"+current_time+"_slc6a4.txt"

#Save some info to log file
log_file = LOG_DIR + "/" + current_time + "_log.txt"
with open(log_file, "w+") as log:
	log.write("Log created on: " + current_time + "\n")
	
	log.write("\n## GENERAL CONFIGURATION ##\n")
	log.write("Population: " + POP + "\n")
	log.write("Additional SNPs: " + ','.join(SNPS) + "\n")
	log.write("Web app folder: " + APP_DIR + "\n")

	log.write("\n## RESOURCES FILES ##\n")
	log.write("CYP2D6 metabolyzer phenos: " + cyp2d6_pheno_file + "\n")
	log.write("CYP2C19 metabolyzer phenos: " + cyp2c19_pheno_file + "\n")
	log.write("CYP2D6 diplotype AF: " + cyp2d6_af_table + "\n")
	log.write("AIF file: " + aif_file + "\n")
	
	log.write("\n## INPUT FILES ##\n")
	log.write("Negative control sample id: " + args.ntc_id + "\n")
	log.write("Sample information file: " + args.sample_data + "\n")
	log.write("Allele typer file: " + args.cyp_alleles + "\n")
	log.write("TaqMan genotyper file: " + args.genos + "\n")
	log.write("SLC6A4 alleles: " + args.slc6a4 + "\n")
	
	log.write("CYP2D6 diplotype not found: " + ';'.join(missed_diplotypes) + "\n")
