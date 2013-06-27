#!/usr/bin/perl -w

use Bio::SeqIO;

#modify the inputfile filename as needed
$inputfile = "SampleFile.txt";

# you shouldn't have to modify anything below this line

#print "hello world\n";

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

# read in the data file
# requires "use Bio::SeqIO" above
#open INFILE, 'SampleFile.txt' or die $! ;
#@lines=<INFILE>;
#close(INFILE);
$seqio_obj = Bio::SeqIO->new(-file => $inputfile, -format => "fasta" );

$seq_obj = $seqio_obj->next_seq;  # for testing 

#print $seq_obj->display_id."\n";
#print $seq_obj->desc."\n";
#print $seq_obj->seq."\n";
print substr $seq_obj->seq,0,4;
exit;

while ($seq_obj = $seqio_obj->next_seq){   
    # print the sequence   
    print $seq_obj->display_id,"\n";
    #print substr($seq_obj->seq
}

