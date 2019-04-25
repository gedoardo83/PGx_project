#This script process results from genotyping experiment,
#prepare tables used in PGx report generations and update online data
#by Edoardo Giacopuzzi, 2019

#Location of the PGx folder. All paths refer to this folder
PGXFOLDER=/data1/PGx_project

#Command line arguments are as follows:
#Files from OpenArray
#1. CSV file results exported from AlleleTyper
#2. Genotype matrix exported from OpenArray software
#Tab-separated files with one line header prepared by user:
#3. patients IDs (col1) and L/S genotype (col2) for SLC6A4 gene
#4. patients IDs (col1) and first drugs used (col2)

ALLELETYPER=$1
GENOMATRIX=$2
SLC6A4=$3
FIRSTDRUG=$4

#Additional resource files used. Do not modify these files!
#ASSAYID: assayID (col1), gene_name (col2) and rsID (col3)
#METATABLE: cyp gene (col1), star allele (col2) and metabolyzer pheno (col3)

ASSAYID=resources/AssayID_to_rsID.csv
METATABLE=resources/Metabolizer_table.csv
HEADER=resources/Header.txt

#Output for each single run is generated automatically in sample_tables folder

CURRENTDATE=`date +"%Y%m%d-%H%M"`
echo $CURRENTDATE
OUTPUT=sample_tables/${CURRENTDATE}_samples_genos.csv

#Step 1.Take input files and run perl script to prepare samples_genos table
perl $PGXFOLDER/Prepare_data.pl $ALLELETYPER $GENOMATRIX $SLC6A4 $FIRSTDRUG $PGXFOLDER/$ASSAYID $PGXFOLDER/$METATABLE $PGXFOLDER/$OUTPUT

#Archive original input files into input_data folder
filename=$(basename $ALLELETYPER)
mv $ALLELETYPER $PGXFOLDER/input_data/${CURRENTDATE}_$filename
filename=$(basename $GENOMATRIX)
mv $GENOMATRIX $PGXFOLDER/input_data/${CURRENTDATE}_$filename
filename=$(basename $SLC6A4)
mv $SLC6A4 $PGXFOLDER/input_data/${CURRENTDATE}_$filename
filename=$(basename $FIRSTDRUG)
mv $FIRSTDRUG $PGXFOLDER/input_data/${CURRENTDATE}_$filename

#Merge currently generated file with previous dataset
#Old dataset is saved ALL_samples_genos.old
if [ -f "$PGXFOLDER/sample_tables/ALL_samples_genos.csv" ]; then
  tail -n+2 $PGXFOLDER/sample_tables/ALL_samples_genos.csv > $PGXFOLDER/TEMP1
  mv $PGXFOLDER/sample_tables/ALL_samples_genos.csv $PGXFOLDER/sample_tables/ALL_samples_genos.old
fi

tail -n+2 $PGXFOLDER/$OUTPUT > $PGXFOLDER/TEMP2

if [ -f "$PGXFOLDER/TEMP1" ]; then
  cat $HEADER $PGXFOLDER/TEMP1 $PGXFOLDER/TEMP2 > $PGXFOLDER/sample_tables/ALL_samples_genos.csv
else
  cat $HEADER $PGXFOLDER/TEMP2 > $PGXFOLDER/sample_tables/ALL_samples_genos.csv
fi

rm $PGXFOLDER/TEMP*

SAMPLECOUNT=`cat $PGXFOLDER/sample_tables/ALL_samples_genos.csv | wc -l`

echo "Date and time:" > $PGXFOLDER/sample_tables/Current_version.log
echo $CURRENTDATE >> $PGXFOLDER/sample_tables/Current_version.log
echo "Number of samples:" >> $PGXFOLDER/sample_tables/Current_version.log
echo $(($SAMPLECOUNT-1)) >> $PGXFOLDER/sample_tables/Current_version.log

#Check if sample information are provided in samples_info.csv file
perl $PGXFOLDER/Check_samples.pl $PGXFOLDER/sample_tables/ALL_samples_genos.csv $PGXFOLDER/samples_info.csv

#Update files into web_app folder
if [ -f "$PGXFOLDER/web_app/samples_genos.csv" ]; then
  mv $PGXFOLDER/web_app/samples_genos.csv $PGXFOLDER/web_app/samples_genos.old
fi
cp $PGXFOLDER/sample_tables/ALL_samples_genos.csv $PGXFOLDER/web_app/samples_genos.csv

if [ -f "$PGXFOLDER/web_app/samples_info.csv" ]; then
  mv $PGXFOLDER/web_app/samples_info.csv $PGXFOLDER/web_app/samples_info.old
fi
cp $PGXFOLDER/samples_info.csv $PGXFOLDER/web_app/samples_info.csv

SAMPLECOUNT=`cat $PGXFOLDER/web_app/samples_genos.csv | wc -l`
echo "Date and time:" > $PGXFOLDER/web_app/Data_version.log
echo $CURRENTDATE >> $PGXFOLDER/web_app/Data_version.log
echo "Number of samples:" >> $PGXFOLDER/web_app/Data_version.log
echo $(($SAMPLECOUNT-1)) >> $PGXFOLDER/web_app/Data_version.log

#Upload app files
cd $PGXFOLDER/web_app
Rscript --vanilla UpdateApp.R

cd ..