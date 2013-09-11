#!/usr/bin/perl -w

# find key sequences in an input file.
# a keyword for incomplete code is TODO

use Bio::SeqIO;

#modify the inputfile filename as needed
$inputfile = "SampleFile.txt";

$offset=91;  # offset between reference sequence and reads

$fudge_factor = 2;  #the number of insertions or deletions allowed before the sequence of interest


##########################################################
# you shouldn't have to modify anything below this line
##########################################################

$debug=0; #set to 1 for some debugging output

$nSilent=20; #number of silent copies, including the reference

$nRegions=10;

# define the regions
@region_min=(0, #{{{
    159,
    187,
    223,
    253,
    273,
    289,
    315,
    347,
    394,
    474);#}}}

@region_length=(0, #{{{ 
    181-159+1, #make the computer do the subtraction
    210-187+1,
    239-223+4,
    273-253+1,    
    290-273+1,
    305-289+1,
    332-315+1,
    365-347+1,
    447-394+25,
    494-474+1);#}}}

@index = ("ref", #{{{
    "1c1",
    "1c2",
    "1c3",
    "1c4",
    "1c5",
    "2c1",
    "2c2",
    "2c3",
    "2c4",
    "2c5",
    "2c6",
    "3c1",
    "3c2",
    "3c3",
    "6c1",
    "6c2",
    "6c3",
    "7c1",
    "uss"); #}}}

@R1=("CGTCACCGAGTATTACCTGAATC", #{{{
    "CGTTGCCGGGTATTGCCCGAATC",
    "CGTTACCGGGTATTGCCCGAATC",
    "CGTTGCCGGGTATTGCCCGAATC",
    "AGTTGCCGGGTATTGCCTGAATC",
    "CGTTACCGAGTATTACCTGAATC",
    "CGTCACCGAATATTACCCGAATA",
    "CGTTACCGAGTATTGCCCGAATC",
    "CGTTACCGAGTATTACCTGAATC",
    "GGTTGCCGGGTATTGCCTGAATC",
    "CGTTGCCGGGTATTACCTGAATC",
    "CGTTACCGGGTATTGCCCGAATC",
    "CGTTGCCGGGTATTACCTGAATC",
    "CCTCACCGAATATTACCTGAAA",
    "CGTCACCGAGTATTACCTGAATC",
    "CGTTACCGGGTATTACCTGAATA",
    "CGTTGCCGGGTATTACCCGAATC",
    "CGTCACCGAGTATTACCCGAATA",
    "CGTCACCGAGTATTACCTGAATC",
    "CGTCACCGAGTATTACCTGAATC"); #}}}

@R2=("ATATGGCCGAAAGACAACACTTCT", #{{{
    "ACATGGCCGAAAGACAACGGTGAT",
    "AAATGGCCGAAAGACAACACTTCT",
    "AAATGGCCGGAAAACAACACTTCT",
    "GAATGGCCGGAAGACAACACTTCT",
    "GAATGGCCCAAAGACAACGGCTCT",
    "GAATGGCCCGCCGACAACGGCGCT",
    "GAATGGCCGAAAGACAACGACAAG",
    "ATATGGCCGGAAAACAACGACAAG",
    "GAATGGCCGGAAGACAACACTTCT",
    "GAATGGCCGGAAAACAACGACAAG",
    "ATATGGCCGAAAGACAACACTTCT",
    "GAATGGCCGGAAAACAACACTTCT",
    "GAATGGCCCAAAGACAACGGCTCT",
    "ATATGGCCGAAAGACAACACTTCT",
    "AAATGGCCCGCCGACAACGGCGCT",
    "AAATGGCCGGAAGACAACACTTCT",
    "AAATGGCCCGCAAACAACGGCGCT",
    "ACATGGCCGGAAAACAACACTTCT",
    "ACATGGCCGAAAGACAACGGTGAT") ;#}}}

@R3=("TCCGCTTCAACAATCAA", #{{{
    "TCCCCCGCCGACAAAATCAA",
    "TCCCCCGCCGAAATCAA",
    "TCCCCCCCCCTCCGACATCAA",
    "TCCCCCCCCACCGACATCAA",
    "TCCGCTTCAAAAATCAT",
    "TCCGCTTCAACAATCAA",
    "TCCCCCCCCTCCAACATCAA",
    "TCTTCTTCATCAATCAA",
    "TCCCCCCCACCGACATCAA",
    "TCCGCCTCCGACATCAA",
    "TCCTCCGCCGCCGACATCAA",
    "TCCTCCGACAAAATCAA",
    "TCCGCTTCAACAATCAA",
    "AACCCCACCGACATCAA",
    "TCCCCCGCCACCGACATCAA",
    "TCCCCCCCTCCGACATCAA",
    "TCCGCCTCCGACATCAA",
    "TCCTCCGCCACCGACATCAA",
    "TCCTCCGCCGAAATCAA"); #}}}

@R4=("CAGAAAGTTGAAGTCGCAAAA", #{{{
    "CAGAAAGTTGAAGTCGCAAAA",
    "AAAAGCGTTACGGTCGCAAAA",
    "AAAGAGGTTGAAGTTAAAAAC",
    "AAAAGCGTTACGGTCGCAAAA",
    "AAGGAAGTTAAAGTCGAAAAC",
    "CAGAAAGTTGAAGTCGCAAAA",
    "GAAAGCGTTACGGTCACAAAC",
    "AAGGAAGTTAAAGTCGAAAAC",
    "AAAAGCGTTACGGTCGCAAAA",
    "AAAAGCGTTACGGTCGCAAAA",
    "CAGAAAGTTGAAGTCAACAAC",
    "CAGAAAGTTGAAGTCGCAAAA",
    "CAGAAAGTTGAAGTCACAAAC",
    "GAAAGCGTTACGGTCACAAAC",
    "AAGGAAGTTAAAGTCGAAAAC",
    "CAAAGCGTTACGGTCGCAAAC",
    "GAAAGCGTTACGGTCGAAAAA",
    "CAAAGCGTTACGGTCGCAAAC",
    "AAAAGCGTTACGGTCGCAAAC"); #}}}

@R5=("AGGCGTCGTTACCGCCCA", #{{{
    "AGGCGTCGTTACCGCCGA",
    "AGGCGTCGTTACCGCCCA",
    "CGGCGTCGTTACCGCCAC",
    "AGGCGTCGTCACCGCCGA",
    "CGGCGTCGTCACCGCCCA",
    "AGGCGTCGTTACCGCCCA",
    "CGGCGTCGTTACCGCCAC",
    "CGGCGTCGTCACCGCCAC",
    "AGGCGTCGTCACCGCCGA",
    "AGGCGTCGTTACCGCCGA",
    "CGGCGTCGTTACCGCCAC",
    "AGGCGTCGTTACCGCCCA",
    "CGGCGTCGTTACCGCCCA",
    "CGGCGTCGTTACCGCCAA",
    "CGGCGTCGTCACCGCCCA",
    "CGGCGTCGTTACCGCCGA",
    "AGGCGTCGTTACCGCCAA",
    "CGGCGTCGTTACCGCCGA",
    "CGGCGTCGTTACCGCCGA"); #}}}

@R6=("CAAATGGCTTCAACCGG", #{{{
    "GAAATGAAACCAAGCGG",
    "CAAATGAATCCAAGCGG",
    "ACAATGCTTTCAAGCGG",
    "GAAATGGCTTCAACCGG",
    "CAAATGGCTTCAAGCAA",
    "CAAATGGCTTCAACCGG",
    "ACAATGCTTTCAAGCGG",
    "ACAATGAATTCAAGCAA",
    "GAAATGGCTTCAACCGG",
    "GAAATGGCTTCAACCGG",
    "ACAATGGCTTCAAGCAA",
    "CAAATGGCTTCAACCGG",
    "CAAATGGCTTCAACCGG",
    "AAAATGCTTTCAAGCGG",
    "CAAATGGCTTCAACCGG",
    "GAAATGAAATCAGACGG",
    "AAAATGCTTTCAAGCGG",
    "GAAATGAAATCAGACGG",
    "GAAATGGCTTCAAGCGG"); #}}}

@R7=("AGAAATCCAAGACAAAAA", #{{{
    "AGAAATCAAAGGCAAAAA",
    "TGAAATCAAAGACAAAAA",
    "TGAAATCAAAGGCAAAAA",
    "TGAAATCAAAGGCAAAAA",
    "AGAAATCAAAGACAAAAA",
    "AGAAATCCAAGACAAAAA",
    "TGAAATCAAAGGCAAAAA",
    "AGAAATCAAAGACAAAAG",
    "TGAAATCAAAGGCAAAAA",
    "AGAAATCCAAGGCAAAAA",
    "AGAAATCCAAGGCAAAAG",
    "AGAAATCCAAGGCAAAAA",
    "AGAAATCAAAGGCAAAAA",
    "AGAAATCCAAGGCAAAAG",
    "TGAAATCAAAGGCAAAAA",
    "AGAAATCAAAGGCAAAAA",
    "AGAAATCCAAGGCAAAAA",
    "AGAAATCCAAGGCAAAAG",
    "AGAAATCAAAGGCAAAAG"); #}}}

@R8=("CCAAGCGTCAAGACGGTTC", #{{{
    "CCAAGCGTGAAGACGGTTC",
    "CCAAGCGTGAAAACGGTTC",
    "CCAGGCGTGAAAACGGTTC",
    "CCAAGCGTCAAGACGGTC",
    "CCAAGCGTGAAAACGGTTC",
    "CCAAGCGTCAAGACGGTTC",
    "CCAAGCGTCAAGCCGGTTC",
    "CCAAGCGTGAAAACGGTTC",
    "CCAAGCGTCAAGACGGTTC",
    "CCAAGCGTGAAAACGGTTC",
    "CCAAGCGTGAAAACGGTC",
    "CCAAGCGTCAAGACGGTTC",
    "CCAGGCGTCAAGACGGTTC",
    "CCAAGCGTGAAGCCGGTTC",
    "CCAAGCGTCAAGACGGTTC",
    "GCAGGCGTGAAAACGGTTC",
    "CCAAGCGTGAAGCCGGTTC",
    "CCAGGCGTGAAGCCGGTTC",
    "CCAGGCGTCAAGACGGTTC"); #}}}

@R9=("ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC", #{{{
    "AAGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA",
    "ACGCGCAACGACGCCAAAGCCGACGCCAAAGACGACACCGTCACCGCCATCGAA",
    "ACGCGCGCCAAAGCCGACGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC",
    "AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAAAGACACCGACAAAATCAAC",
    "AAGCGCACCGAAGCCAACGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACACCGCCGGCACCAAAATCGAC",
    "ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC",
    "CAGCGCGCCAAAGCCGACGACGCCGTCACCGCCGACGCCAACAACGCCATCGAC",
    "AAGCGCGCCAACGTTGCCGCCGCCAACGACGACGACGTTACCGACGACAAAAACAACAACGGCATCGAC",
    "AAGCGCGCCGACAACAACGGCAACATTACCGCCGACAACGGCAACGCCATCGAA",
    "AAGCGCACCGAAGCCAACGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACACCGCCGGCACCAAAATCGAC",
    "ACGCGCGACGACAAAGCCAAAGACGACGTCAAAGCCGACGGCACCGCCGGCACCAAAATCGAC",
    "ACGCGCAACGACGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA",
    "AAGCGCGACGACGCCGCCGCCAAAGACGACACCGTCACCGCCGACGCCACCGGCAACGACGGCAAAATCGAC",
    "AAGCGCACCGAAGCCAACGCCGACGCCGCCGGCAAAGACACCACCAACGGCATCAAC",
    "AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAACAACGGCATCAAC",
    "AAGCGCGACGCCAACAACGCCAACAACGACGCCGTCACCGACGACACCACCGGCAACGGCAACGAAAAAATCGAA",
    "AAGCGCAACGACGCCGCCAACGACGACGTTACCGACGACGCCGGCACCGACAACGGCGGCAAAGGCAAAATCGAC",
    "ACGCGCGCCAAAGCCAAAGACGCCGACGACGTTACCGACGACGCCGGCACCCACAACGGCGGCAAAGGCAAAATCGAC",
    "ACGCGCAACGACGCCGCCGACAACGACGACGTCGCCAAAGACGACGCCGCCGGCAACGCCATCGAA"); #}}}

@R10=("CGATGAATCATCTGCCACCTA", #{{{
    "CGATGAATCATCTGCCACCTA",
    "CGATGAATCATCTGCCACCTA",
    "CGATGAATCATCTGCCGTT",
    "CGATAAATCATCTGCCGTT",
    "CGATGAATCATCGTTGCCGG",
    "CGATAAACATGATGCCAAATG", 
    "TGATACGTCATCTGCCAAA",
    "CGATAAATCATCTGCCACCTA",
    "CGATGAATCATCTGCCACCTA",
    "CGATGAATCATCGTTGCCGG",
    "CGATGAATCATCGTTGCCGG",
    "CGATGAATCATCTGCCACCTA",
    "CGATAAATCAACTGCCGTT",
    "CGACCCGTTCTCTGCTAGC",
    "CGATAAACATGATGCCAAA",
    "CGATGAATCATCTGCCGTTTA",
    "CGATAAATCAACTGCCGTT",
    "CGATAAATCAACTGCCAAA",
    "CGATGAACCAACTGCCACCTA"); #}}}

# simplify the region search pointer {{{
$R[1]=\@R1;
$R[2]=\@R2;
$R[3]=\@R3;
$R[4]=\@R4;
$R[5]=\@R5;
$R[6]=\@R6;
$R[7]=\@R7;
$R[8]=\@R8;
$R[9]=\@R9;
$R[10]=\@R10;#}}}

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

if ($debug>=2){ # for testing 
  for ($i=0;$i<10;$i++){
    $seq_obj = $seqio_obj->next_seq;
  } 
  print $seq_obj->display_id."\n";
  print $seq_obj->desc."\n";
  print $seq_obj->seq."\n";
  print "number of basepairs ".length $seq_obj->seq;
  print "\n";

  for ($r=1; $r<=$nRegions; $r++){ #start later for debugging
    for ($i=0; $i<$nSilent; $i++){
      if ($debug>=2) {
	print "\$r=$r, \$i=$i\n"
      }
      my $substr = substr $seq_obj->seq, $region_min[$r]-$offset-$fudge_factor, $region_length[$r]+2*$fudge_factor;
      my $search = $R[$r][$i];
      if ($debug>=2) {
	print "\$substr=".$substr."\n";
	print "\$search=".$search."\n";
      }
      if ( $substr =~ /$search/ ){
	print "FOUND ";
	print $index[$i];
#print " at position ";
#print pos($substr);
	print "\n";
      } 
    } 
  }


  exit;
}

#header
print "counter,ID,region1,region2,region3,region4,region5,region6,region7,region8,region9,region10\n";

$counter=1;
while ($seq_obj = $seqio_obj->next_seq) {   
    print $counter,",";
    $counter+=1;
    #if ($counter==4) {exit;}
    # print the display_id   
    print $seq_obj->display_id,",";

    #check sequence is long enough to run the tests
    if (length($seq_obj->seq) < $region_min[$nRegions]+$region_length[$nRegions]){ #+2*$fudge_factor){
        print "error: sequence too short\n";
        next;
    }

    for ($r=1; $r<=$nRegions; $r++){ #start later for debugging
        for ($i=0; $i<$nSilent; $i++){
            if ($debug>=1) {
                print "\$r=$r, \$i=$i\n"
            }
            my $substr = substr $seq_obj->seq, $region_min[$r]-$offset-$fudge_factor, $region_length[$r]+2*$fudge_factor;
            my $search = $R[$r][$i];
            if ($debug>=1) {
                print "\$substr=".$substr."\n";
                print "\$search=".$search."\n";
            }
            if ( $substr =~ /$search/ ){
                print $index[$i]," ";
                if ($index[$i] eq "ref"){
                    last;  # if it's ref, don't bother checking the others
                }
            } 
        } 
        print ",";
    }
    print "\n";
}
