#!/usr/bin/perl -w

use Bio::SeqIO;

#modify the inputfile filename as needed
$inputfile = "SampleFile.txt";

# you shouldn't have to modify anything below this line
$debug=1; #set to 1 for some debugging output


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
if ($debug==2) {print $ref_seq."\n";}


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
if ($debug==2){
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
  for ($i=0;$i<10;$i++){
    $seq_obj = $seqio_obj->next_seq;
  } 
  print $seq_obj->display_id."\n";
  print $seq_obj->desc."\n";
  print $seq_obj->seq."\n";
  print substr $seq_obj->seq,60+60+27,9;
  print "\n";
  print "  ",substr $seq_obj->seq,60+60+29,5;
  print "\n";
  print "    ",substr $seq_obj->seq,60+60+31,1;
  print "\n";
  print "exiting after 1 read\n";
  exit;
}

# do the tests
print "display Id,test base 239,test base 464\n"; #header

$offset=88; #normally 88
$counter=1;
while ($seq_obj = $seqio_obj->next_seq) {   
  print $counter,",";
# print the display_id   
  print $seq_obj->display_id,",";

# test pair 239 (the offset on all reads is 90) 
  $base239 = substr $seq_obj->seq,239-$offset,9;
# is this a useful output, or would true/false be better?
  print $base239;
  print ",";
# test pair 239 (the offset on all reads is 90)
  $base474 = substr $seq_obj->seq,474-$offset,9;
  print $base474;
  print "\n";
  $counter+=1;
}

