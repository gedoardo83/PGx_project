sed -i 's/^M/\n/g' CYP_alleles.csv
tail -n+11 CYP_alleles.csv > CYP_allels.tmp
mv CYP_alleles.tmp CYP_alleles.csv
perl 
