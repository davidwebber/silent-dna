#!/usr/bin/perl -w

use Bio::SeqIO;

#modify the inputfile filename as needed
$inputfile = "SampleFile.txt";

# you shouldn't have to modify anything below this line
$debug=0; #set to 1 for some debugging output


#read in the reference sequence
open INFILE, 'Reference seq.txt' or die $! ;
@lines=<INFILE>;
close(INFILE);
foreach $line (@lines) {
  if ($line =~ /ATGAA/ ){
    chomp $line;
    $ref_seq=$line;
  }
}
if ($debug==1) {print $ref_seq."\n";}


#read in the variants
open INFILE, 'FA1090 MS11 silent copies.txt' or die $! ;
@lines=<INFILE>;
close(INFILE);

$i=0;
$j=0;
foreach $line (@lines) {
  if ( $line =~ /pil/ ){ # if "pil" appears in the line
    #print $line; #debugging
    chomp $line;
    $variant_names[$i++]=$line;
  }
  if ($line =~ /A/ ){
    chomp $line;
    $variant_seq[$j++]=$line;
  }
}
if ($debug==1){
  print @variant_names;
  foreach $line (@variant_seq){
    print $line."\n";
  }
}
if ($i != $j){
  print "uh oh, there's a mismatch between the names of the variants and the sequences\n";
}

# read in the data file
$seqio_obj = Bio::SeqIO->new(-file => $inputfile, -format => "fasta" );

if ($debug==1){ # for testing 
  $seq_obj = $seqio_obj->next_seq; 
  print $seq_obj->display_id."\n";
  print $seq_obj->desc."\n";
  print $seq_obj->seq."\n";
  print "exiting after 1 read\n";
  exit;
}

# do the tests
print "display Id,test 239,test 464\n"; #header

while ($seq_obj = $seqio_obj->next_seq) {   
# print the display_id   
  print $seq_obj->display_id,",";

# test pair 239 (the offset on all reads is 91) 
  $base239 = substr $seq_obj->seq,239-91,1;
# is this a useful output, or would true/false be better?
  print $base239;
  print ",";
# test pair 239 (the offset on all reads is 91)
  $base474 = substr $seq_obj->seq,474-91,1;
  print $base474;
  print "\n";
}

