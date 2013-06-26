#!/usr/bin/perl -w

# note to self: use bio seqIO

#print "hello Ella\n";

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
#print $ref_seq."\n";


#read in the variants
open INFILE, 'FA1090 MS11 silent copies.txt' or die $! ;
@lines=<INFILE>;
close(INFILE);

$i=0;
$j=0;
foreach $line (@lines) {
# if "pil" appears in the line
  if ( $line =~ /pil/ ){
    #print $line;
    chomp $line;
    $variant_names[$i++]=$line;
  }
  if ($line =~ /A/ ){
    chomp $line;
    $variant_seq[$j++]=$line;
  }
}
#print @variant_names;
#print @variant_seq;
if ($i != $j){
  print "uh oh, there's a mismatch between the names of things and the sequences\n";
}


#read in the data file
open INFILE, 'SampleFile.txt' or die $! ;
@lines=<INFILE>;
close(INFILE);



