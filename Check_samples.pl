#Part of the Process_results.sh pipeline
#Check if sample information are present in samples_info.csv file for all samples with genotypes

use Array::Utils qw(:all);

$genofile=$ARGV[0];
$infofile=$ARGV[1];

open(IN, $genofile);
$header=<IN>;
while($row=<IN>) {
  chomp($row);
  @line=split("\t", $row);
  push(@genosamples, $line[0]);
}

open(IN, $infofile);
$header=<IN>;
while($row=<IN>) {
  chomp($row);
  @line=split("\t", $row);
  push(@infosamples, $line[0]);
}

@diff = array_minus(@genosamples, @infosamples);

if (scalar @diff > 0) {
  $message = << "END_MESSAGE";
  WARNING!! scalar(@diff) SAMPLES ARE MISSING FROM sample_info.csv
  Following samples have genotypes but are not found:
  join("\n",@diff)
END_MESSAGE
  print $message;
}