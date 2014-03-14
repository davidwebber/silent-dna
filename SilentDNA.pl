#!/usr/bin/perl -w

# find key sequences in an input file.
# a keyword for incomplete code is TODO

# To run: ./SilentDNA.pl inputfile.fna [ first_read ] [ last_read ] > outputfile.csv

# sequences are labeled according to specific silent copy names.
# if a sequence is tagged but has no specific silent copy, the keyword is "var"

# DONE   tolerate short reads
# DELETE implement short regions in region 9 
# DONE   implement a "verdict" column
# TODO   try to index region 10 from the right, since R9 length varies
# TODO   count the number of empty cells on each row
# TODO   columns for counts of "var (any incl named)", "blank", and "ref"

#use Text::LevenshteinXS qw(distance);
use Text::WagnerFischer qw(distance);

use Bio::SeqIO;
use warnings;
no warnings ('substr');
#no warnings ('uninitialized', 'substr');
require 'SequenceData1.pl';

$debug=0; #set to 1 or 2 for some debugging output

$inputfile = $ARGV[0];
if (not defined $inputfile) {
    die "need input filename";
}

$first_read = $ARGV[1];
if (not defined $first_read){
  $first_read = 1;
}

$last_read = $ARGV[2];
if (not defined $last_read){
  $last_read = 1e25; # a big number
}

$fudge_factor = 3;  #the number of insertions or deletions allowed before the sequence of interest

$boundary_search_width=10; #the allowed fudge from where we search for the boundary and where we expect it, except region 10.


##########################################################
# you shouldn't have to modify anything below this line
##########################################################

if ($debug>0){
  if ("Hello World" =~ /Hello/){
    print "dude, pattern matching worked!\n";
  }
}

if ($debug>0){
  print $R[2][0]."\n";
}

# read in the data file
$seqio_obj = Bio::SeqIO->new(-file => $inputfile, -format => "fasta" );

#if ($debug>=2){ # for testing 
#  for ($i=0;$i<10;$i++){
#    $seq_obj = $seqio_obj->next_seq;
#  } 
#  print $seq_obj->display_id."\n";
#  print $seq_obj->desc."\n";
#  print $seq_obj->seq."\n";
#  print "number of basepairs ".length $seq_obj->seq;
#  print "\n";
#
#  for ($r=1; $r<=$nRegions; $r++){ #start later for debugging
#    for ($i=0; $i<$nSilent; $i++){
#      if ($debug>=2) {
#	    print "\$r=$r, \$i=$i\n"
#      }
#      my $substr = substr $seq_obj->seq, $region_min[$r]-$offset-$fudge_factor, $region_length[$r]+2*$fudge_factor;
#      my $search = $R[$r][$i];
#      if ($debug>=2) {
#	    print "\$substr=".$substr."\n";
#	    print "\$search=".$search."\n";
#      }
#      if ( $substr =~ /$search/ ){
#	    print "FOUND ";
#	    print $index[$i];
#         #print " at position ";
#         #print pos($s#ubstr);
#	    print "\n";
#      } 
#    } 
#  }
#  exit;
#}

sub printNvar{
    my @result = @_;
    my $total=0;
    my $i=1;
    for ($i=1; $i<=$nRegions; $i++){
        if ($result[$i] ne " " && $result[$i] ne "ref" && $result[$i] ne "short"){
            $total++;
        }
    }
    print $total.",";
}

sub printNblank{
    my @result = @_;
    my $total=0;
    my $i=1;
    for ($i=1; $i<=$nRegions; $i++){
        if ($result[$i] eq " "){
            $total++;
        }
    }
    print $total.",";
}

sub getNblank{
    my @result = @_;
    my $total=0;
    my $i=1;
    for ($i=1; $i<=$nRegions; $i++){
        if ($result[$i] eq " "){
            $total++;
        }
    }
    return $total;
}



sub getVerdict{
    my @result = @_;
    my $resultString = join(",",@result);
    #print $resultString."\n";
    my $nRef =()= $resultString =~ /ref/gi;
    if ($nRef==10){
        #if all are ref, print ref
        return "ref";
    } elsif ($resultString =~ /var/){
        #if it's var, just print var
        return "var";
    } elsif ($resultString =~ /N\/A/){
        #if it's var, just print var
        return "trunc";
    } elsif (getNblank(@result) != 0){
        return ""; #blank
    } else {
        # map each silent copy name to binary and do the bitwise and;
        $totalBitString = $bitmap{'ref'}; # start with all ones
        for (my $i=1; $i<=$nRegions; $i++){
            $subBitString = 0;
            @tokens = split(' ',$result[$i]);
            #TODO but what if $results[$i]=' '?  What is the token then?
            $ntokens = @tokens;
            if ($ntokens!=0){  # OR each bitmap together
            	foreach $token (@tokens) {
                    if (length($token)==3) {
                        $subBitString = $subBitString | $bitmap{$token};
                        #print $token."\n";
                    } 
                    if ($token =~ /short/) {
                      $subBitString = $bitmap{'ref'};# set to all ones 
                    }
                    if ($token =~ /false/) {
                      $subBitString = $bitmap{'ref'};# set to all ones 
                    }
                    if ($token =~ /indel/) {
                      $subBitString = $bitmap{'ref'};# set to all ones 
                    }
                    if ($token =~ /mismatch/) {
                      $subBitString = $bitmap{'ref'};# set to all ones 
                    }
            	}
            } else { #if the cell was blank,
                #$subBitString = $bitmap{'ref'};# set to all ones	
                #$subBitString = $bitmap{'var'};# set to all zeros  
                $subBitString = $bitmap{'should not get here'};# set to all zeros  
            }
            #printf ("\n%20b",$subBitString);
            #print " ".$result[$i]."\n";
            #print "\n".$tokens[1];
            $totalBitString = $totalBitString & $subBitString;
        }
        #printf ("%20b,",$totalBitString);
        if (exists $reversebitmap{$totalBitString}){
            return $reversebitmap{$totalBitString};
        } else {
            return "var";
        }
    }
}

sub printPminus{
    my @result = @_;
    if ( ($result[1] =~ /3c2/) || ($result[3] =~ /2c4/) || ($result[3] =~ /1c3/) || ($result[3] =~ /6c2/) ) {
        print "P-,";
    } else {
        print "-,";
    }
}

sub getAltVerdict{
    # if a variant appears next to a named silent copy, treat that region as ref to allow the silent copy to be the verdict
    my @result = @_;
    #print "\n";
    #for ($r=1;$r<=$nRegions;$r++){
    #        print $result[$r]." ,";
    #};
    #print "\n";
    for (my $i=2; $i<$nRegions; $i++){ # search inners
        if ($result[$i] =~ /v/ && ($result[$i-1] =~ /c/ || $result[$i+1] =~ /c/) ){
            $result[$i] = "ref";
        }       
    }
    if ($result[2] =~ /v/ && $result[1] =~ /c/){ #left side
        $result[2] = "ref"
    }
    if ($result[$nRegions] =~ /v/ && $result[$nRegions-1] =~ /c/){ #right side
        $result[$nRegions] = "ref"
    }
    #print "\n";
    #for ($r=1;$r<=$nRegions;$r++){
    #        print $result[$r]." ,";
    #};
    #print "\n";
    return getVerdict(@result);
}


sub min_distance{
    my ($seq, $ref) = @_;
    my $n=0;
    my $m=0;
    my $length = length $seq;
    my $dist_min=100;

    if ($length < length($ref)){ #short read
        return -1;
    }
 
    my $last_dist;
    my $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
    #print $dist."\n";

    # I think this algorithm could be improved
    for ($n=1; $n<10; $n++){
        $last_dist = $dist;
        $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
        #print "$n $m $dist"."\n";
        if ($dist <= $last_dist){
            $last_dist=$dist;
        } else {
            $n-=1;
            $dist = $last_dist;
            last; #break;
        }
    }

    for ($m=1; $m<10; $m++){
        $last_dist = $dist;
        $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
        #print "$n $m $dist"."\n";
        if ($dist <= $last_dist){
            $last_dist=$dist;
        } else {
            $m-=1;
            $dist = $last_dist;
            last; #break;
        }
    }

    #print "minimum distance is $dist \n";
    return $dist;

}



#header
print "counter,ID,region1,region2,region3,region4,region5,region6,region7,region8,region9,region10,nVar,nEmpty,verdict,altVerdict,P_minus,basepairs\n";
my @result;

$counter=1;

while ($counter < $first_read){
  $seq_obj = $seqio_obj->next_seq or die "file is shorter than first requested sequence";
  $counter+=1;
}

while ( ($seq_obj = $seqio_obj->next_seq) && ($counter<=$last_read) ) {   

    $length_shift = 0; # accumulated shift in the position from long or short regions
    $this_offset = $offset; # this offset may change eg if it's a short read
    for ($r=0; $r<=$nRegions; $r++){
        $result[$r]=" ";
    }

  
    print $counter,",";
    $counter+=1;
    #if ($counter==4) {exit;}
    # print the display_id   
    print $seq_obj->display_id,",";

    $rStart=1;

    #check if sequence is a short read.  We can recover if it's 310 to 340 basepairs (it's probably 6c1)
    if (310 < length($seq_obj->seq) && length($seq_obj->seq) < 340){
        $this_offset=271;
        print "short,short,short,short,short,";
        $rStart=6;
        for ($r=1; $r<6; $r++){
            $result[$r]="short";
        }
    }

    #check sequence is long enough to run the tests
    #if (1==2 && length($seq_obj->seq) < $region_min[$nRegions]+$region_length[$nRegions]){ #+2*$fudge_factor){
    #    #print "error: sequence too short\n";
    #    $this_offset=271;
    #    print "short,short,short,short,short,";
    #    $rStart=6;
    #    for ($r=1; $r<6; $r++){
    #        $result[$r]="short";
    #    }
    #    #TODO continue testing the latter part of the sequence
    #    #next; #go to next sequence
    #}


    for ($r=$rStart; $r<=$nRegions; $r++){ 
        my $substr;
        if ( ($r==2) || ($r==3) || ($r==4) || ($r==7) || ($r==8) || ($r==9) || ($r==10) )#index region 10 from the edge of 9 and 10
        {
            #find the match position
            if ($seq_obj->seq =~ /$bound[$r]/) {
                if ( $debug>=1 ){
                    print "region ".($r-1)." to ".$r." boundary found between @- and @+. ";
                    print "region ".$r." begins at ".($region_min[$r]-$this_offset-$fudge_factor+$length_shift).". delta=";
                }
                #use for small adjustments only
                if (abs("@+"-$region_min[$r]) < $boundary_search_width && $r!=10){
                   $length_shift = "@+" - ($region_min[$r]-$this_offset-1); #assignment, not modification;
                }
                #region 10 boundary allowed to scan over the whole sequence
                elsif ( $r==10 ){
                   $length_shift = "@+" - ($region_min[$r]-$this_offset-1); #assignment, not modification;
                }
                if ( $debug>=1) {
                    print $length_shift."\n";
                }
            }
        } 
        $substr = substr $seq_obj->seq, $region_min[$r]-$this_offset-$fudge_factor+$length_shift, $region_length[$r]+2*$fudge_factor;

        if (defined $substr){

            while (my ($key_seq,$silent_copy) = each $R{$r}){

                if ($debug>=2) {
                    print "\$r=$r, $silent_copy\n"
                }
                # substr source start length
 
                my $search = $key_seq;
                if ($debug>=2) {
                    print "\$substr=".$substr."\n";
                    print "\$search=".$search."\n";
                }

                if ( $substr =~ /$search/ ){
                    print $silent_copy," ";
                    $found_end = "@+" ; #MAGIC why do I need quotes here?
                    #$found_begin = "@-" ;
                    $expected_end = ($fudge_factor + length($search));
                    if ( ($r==1) && ($debug>=1) ){print "region 1 found between @- and @+.  expected between $fudge_factor and ".($fudge_factor + length($search))."\n";}
                    #if ( ($r==1) && ($debug>=1) ){print "region 1 found between $found_begin and $found_end.  expected between $fudge_factor and $expected_end. Diff=".($found_end-$expected_end)."\n";}
                    $length_shift += ($found_end - $expected_end);  #shift the search region if the previous region is found outside where it was expected.
                    $result[$r]=$silent_copy;
                    $length_shift += length($search) - $ref_length[$r];
                    if ($debug>=1) {
                        print "\$substr=".$substr."\n";
                        print "\$search=".$search."\n";
                        print @-,"-",@+," ";
                        print $length_shift."\n";
                    }
                        #if ($silent_copy eq "ref"){
                        #   last;  # if it's ref, don't bother checking the others
                        #}
                }        
            } 
            if ($result[$r] eq " "){ # " " is the initial value
                $min_dist = min_distance($substr,$R{$r}{'ref'});
                if ($min_dist == 1) {
                    print "mismatch";
                    $result[$r]="mismatch";
                } else {
                    print $min_dist;
                    print " ";
                    print $substr; # if no match found, print the search sequence to figure out why
                }
            }
        } else { # cell does not exist due to short read
            print "N/A";
            $result[$r]="N/A";
        }
        print ",";
    }

    &printNvar(@result); # print the number of variants
    &printNblank(@result); # print the number of blanks
    $verdict = &getVerdict(@result); # get the verdict
    print $verdict.",";
    $verdict = &getAltVerdict(@result);
    print $verdict.",";
    &printPminus(@result);

    print length($seq_obj->seq); #debug
    
    print "\n"; # end of line

    
    if ($debug>=1){
        print "                 ,";
        for ($r=1;$r<=$nRegions;$r++){
            print $result[$r]." ,";
        }
        print "\n";
    }

}
