#Script to prepare table used in PGx report generations
#Command line arguments are as follows:
#Files from OpenArray
#1. CSV file results exported from AlleleTyper
#2. Genotype matrix exported from OpenArray software
#Tab-separated files with one line header:
#3. patients IDs (col1) and L/S genotype (col2) for SLC6A4 gene
#4. patients IDs (col1) and first drugs used (col2)
#5. assayID (col1), gene_name (col2) and rsID (col3)
#6. cyp gene (col1), star allele (col2) and metabolyzer pheno (col3)
#7. Output file
#
#The script performs the following tasks:
#1. Translate CYP SNPs to metabolyzer phenotype
#2. Translate assayID to rsID
#3. Provide final table with metabolyzer phenos, the SNPs needed for PGx report and L/S genotype

########################################################
###  Define relevant SNPs besides CYP2D6 and CYP2C19 ###
########################################################

@mysnps = ("rs3211371", "rs56337013", "rs35694136"); # ("rs489693","rs4713916","rs7997012","rs6295");

#####################################
###  READ COMMAND-LINE ARGUMENTS  ###
#####################################
#Output file
$outfile = $ARGV[6];

#Table of results from AlleleTyper
$inputfile = $ARGV[0];

#Genotype matrix from OpenArray software
$genofile = $ARGV[1];

#Table with L/S genotypes for SLC6A4
$LSfile = $ARGV[2];

#Table with first drug used for the patients
$drugfile = $ARGV[3];

#Assay id to rsID conversion table (tab-separated file with 3 cols: assayID,gene,rsID)
$assayIDfile = $ARGV[4];

#Assay id to rsID conversion table (tab-separated file with 3 cols: assayID,gene,rsID)
$metabfile = $ARGV[5];

#open output file for writing
open(OUT, ">>$outfile");

print	"##Arguments as interpreted##\n".
		"AlleleTyper file:\t$inputfile\n".
		"Genotype matrix:\t$genofile\n".
		"SLC6A4 L/S:\t$LSfile\n".
		"First drug used:\t$drugfile\n".
		"Assay ID to rsID:\t$assayIDfile\n".
		"metabolyzer phenos:\t$metabfile\n".
		"Output file:\t$outfile\n";

##########################################
###  READ assayID, drug and L/S FILES  ###
##########################################

#Read Assay ID conversion table
print "INFO\tReading Assay ID file...\n";

open(IN, $assayIDfile);
$header = <IN>;

while($row = <IN>) {
  chomp($row);
  @line = split("\t", $row);
  $assayID{$line[0]}=$line[2];
}
close(IN);

#Read SLC6A4 L/S file
print "INFO\tReading SLC6A4 L/S file...\n";

open(IN, $LSfile);
$header = <IN>;

while($row = <IN>) {
  chomp($row);
  @line = split("\t", $row);
  $LSgeno{$line[0]}=$line[1];
}
close(IN);

#Read drug file
print "INFO\tReading drug file...\n";

open(IN, $drugfile);
$header = <IN>;

while($row = <IN>) {
  chomp($row);
  @line = split("\t", $row);
  $drug{$line[0]}=$line[1];
}
close(IN);

#Read metabolyzer phenotype file
print "INFO\tReading metabolyzer file...\n";

open(IN, $metabfile);
$header = <IN>;

while($row = <IN>) {
  chomp($row);
  @line = split("\t", $row);
  $metabpheno{$line[0]}{$line[1]}=$line[2];
}
close(IN);

################################
### PROCESS GENOTYPE MATRIX  ###
################################

print "INFO\tReading genotype matrix file...\n";
open(IN, $genofile);
while($row = <IN>) {
	next if ($row=~ /^#/); #skip comment lines
	next if ($row=~ /^\s*$/); #skip empty lines
	chomp($row);
	$row =~ s/"//g;
	@line = split(",", $row);
	if ($line[0] =~ /^Sample/) { #Process header line
			$i=1;
		foreach(@line[1..$#line]) {
			$cols_index{$assayID{$_}}=$i;
			$i++; 
		}
	} else { #Process genotype lines
		foreach(@mysnps) {
			if ($line[$cols_index{$_}] =~ /[ACTG]\/[ACTG]/) {
				push(@{$output{$line[0]}{genopart}}, $line[$cols_index{$_}]);
			} else {
				push(@{$output{$line[0]}{genopart}}, "NA");	
			} 
		}
	}
}
close(IN);


##################################################################
###  Convert CYP2D6/CYP2C19 alleles to Metabolizer phenotypes  ###
##################################################################
print "INFO\tReading AlleleTyper results file...\n";

#Eventually change ^M to standard \n end-line charcter
$command = 'sed -i \'s/\r/\n/g\' '.$inputfile;
print "$command\n";
system($command);

#Read AlleleTyper results
open(IN, $inputfile);

#Convert CYP alleles to metabolizer phenotypes
while($row = <IN>) {
	next if ($. <= 10); # Skip comment lines
	chomp($row);
	@line = split("\t", $row);
	if ($line[0] =~ /^sample/) { #Process header line
		$i=1;
		foreach(@line[1..$#line]) {
			$cols_index{$_}=$i;
			$i++; 
		}
	} else {
		if ($line[0] eq "Empty" || $line[0] eq "NTC") {next}
		#Process for CYP2D6
			undef @metab;
			@mystring = split(",", $line[$cols_index{CYP2D6}]);
			$mystring[0] =~ /(\*[1-9A-Zx]+)\/(\*[1-9A-Zx]+)/;
			$allele1 = $1;
			$allele2 = $2;
			
			if ($allele1 =~ /(\*[1-9A-Z]+)x([0-9])/) {
				$value = $metabpheno{CYP2D6}{$1};
				push (@metab, ($value) x $2);
			} else {push (@metab, $metabpheno{CYP2D6}{$allele1})}
			if ($allele2 =~ /(\*[1-9A-Z]+)x([0-9])/) {
				$value = $metabpheno{CYP2D6}{$1};
				push (@metab, ($value) x $2);
			} else {push (@metab, $metabpheno{CYP2D6}{$allele2})}
			
			my $UM = grep { $_ eq "UM"  } @metab;
			my $PM = grep { $_ eq "PM" } @metab;
			my $IM = grep { $_ eq "IM" } @metab;
			my $EM = grep { $_ eq "EM" } @metab;

			if ($PM>=2 && $UM==0 && $IM==0 && $EM==0) {
				$output{$line[0]}{CYP2D6}{pheno}="PM";
			} elsif (($PM>=1 && $UM==0 && ($IM>=1 || $EM>=1)) || ($PM==0 && $EM==0 && ($IM==1 || $IM==2)) || ($PM==0 && $IM==1 && $EM==1) || ($PM==0 && $IM==0 && $EM==1)) {
				$output{$line[0]}{CYP2D6}{pheno}="IM";
			} elsif ($EM==2 && $PM==0 && $UM==0) {
				$output{$line[0]}{CYP2D6}{pheno}="EM";
			} elsif (($IM + $EM) > 2) {
				$output{$line[0]}{CYP2D6}{pheno}="UM";
			} else {$output{$line[0]}{CYP2D6}{pheno}="NA"}
			
			$output{$line[0]}{CYP2D6}{alleles}=$mystring[0];

		#Process for CYP2C19
			undef @metab;
			@mystring = split(",", $line[$cols_index{CYP2C19}]);
			$mystring[0] =~ /(\*[1-9A-Zx]+)\/(\*[1-9A-Zx]+)/;
			$allele1 = $1;
			$allele2 = $2;

			if ($allele1 =~ /(\*[1-9A-Z]+)x([0-9])/) {
				$value = $metabpheno{CYP2C19}{$1};
				push (@metab, ($value) x $2);
			} else {push (@metab, $metabpheno{CYP2C19}{$allele1})}
			if ($allele2 =~ /(\*[1-9A-Z]+)x([0-9])/) {
				$value = $metabpheno{CYP2C19}{$1};
				push (@metab, ($value) x $2);
			} else {push (@metab, $metabpheno{CYP2C19}{$allele2})}
			print join("\t", @metab)."\n";

			my $UM = grep { $_ eq "UM"  } @metab;
			my $PM = grep { $_ eq "PM" } @metab;
			my $IM = grep { $_ eq "IM" } @metab;
			my $EM = grep { $_ eq "EM" } @metab;
			if ($PM>=2 && $UM==0 && $IM==0 && $EM==0) {
				$output{$line[0]}{CYP2C19}{pheno}="PM";
			} elsif ((($IM==1 || $IM==2) && $PM <= 1 && $EM==0) || ($EM==1 && $PM <=1)) {
				$output{$line[0]}{CYP2C19}{pheno}="IM";
			} elsif ($EM==1 && $PM==0 && $IM==1) {
				$output{$line[0]}{CYP2C19}{pheno}="EM";
			} elsif ($EM==2 || (($EM + $IM) > 2)) {
				$output{$line[0]}{CYP2C19}{pheno}="UM";
			} else {$output{$line[0]}{CYP2C19}{pheno}="NA"}
			$output{$line[0]}{CYP2C19}{alleles}=$mystring[0];
	}
}
close(IN);

#########################
###  Generate output  ###
#########################

#Generate header
$header="Sample\tCYP2D6_pheno\tCYP2D6_alleles\tCYP2C19_pheno\tCYP2C19_alleles\tSLC6A4\t".join("\t", @mysnps)."\t1st_drug\n";
print OUT $header;

foreach $sample(keys %output) {
	print OUT $sample."\t".
			$output{$sample}{CYP2D6}{pheno}."\t".$output{$sample}{CYP2D6}{alleles}."\t".
			$output{$sample}{CYP2C19}{pheno}."\t".$output{$sample}{CYP2C19}{alleles}."\t".
			$LSgeno{$sample}."\t".join("\t", @{$output{$sample}{genopart}})."\t".$drug{$sample}."\n";
}

close(OUT);
