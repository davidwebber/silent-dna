#!/usr/bin/perl -w

# find key sequences in an input file.
# a keyword for incomplete code is TODO

# TODO tolerate short reads
# TODO implement short regions in region 9 
# TODO implement a "verdict" column
# TODO try to index region 10 from the right, since R9 length varies

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

@index=("ref", #{{{
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

# Region 1 
#$R{1}{'CGTTGCCGGGTATTGCCCGAATC'}='1c1';
#$R{1}{'CGTTACCGGGTATTGCCCGAATC'}='1c2';
$R{1}{'CGTTGCCGGGTATTGCCCGAATC'}='1c1 1c3';
$R{1}{'AGTTGCCGGGTATTGCCTGAATC'}='1c4';
#$R{1}{'CGTTACCGAGTATTACCTGAATC'}='1c5';
$R{1}{'CGTCACCGAATATTACCCGAATA'}='2c1';
$R{1}{'CGTTACCGAGTATTGCCCGAATC'}='2c2';
$R{1}{'CGTTACCGAGTATTACCTGAATC'}='1c5 2c3 var162';
$R{1}{'GGTTGCCGGGTATTGCCTGAATC'}='2c4';
#$R{1}{'CGTTGCCGGGTATTACCTGAATC'}='2c5';
$R{1}{'CGTTACCGGGTATTGCCCGAATC'}='1c2 2c6';
$R{1}{'CGTTGCCGGGTATTACCTGAATC'}='2c5 3c1';
$R{1}{'CCTCACCGAATATTACCTGAAA' }='3c2';
#$R{1}{'CGTCACCGAGTATTACCTGAATC'}='3c3';
$R{1}{'CGTTACCGGGTATTACCTGAATA'}='6c1';
$R{1}{'CGTTGCCGGGTATTACCCGAATC'}='6c2';
$R{1}{'CGTCACCGAGTATTACCCGAATA'}='2c1 6c3';
#$R{1}{'CGTCACCGAGTATTACCTGAATC'}='7c1';
#$R{1}{'CGTCACCGAGTATTACCTGAATC'}='uss';
$R{1}{'CGTCACCGAGTATTACCTGAATC'}='ref'; # uss, 7c1, 3c3
$R{1}{'CGTTGCCGAGTATTACCTGAATC'}='var163';
$R{1}{'CGTCACCGAGTATTACCTGAATA'}='2c1 6c1 6c3 (C181A)';

# Region 2
#$R{2}{'ACATGGCCGAAAGACAACGGTGAT'}='1c1';
$R{2}{'AAATGGCCGAAAGACAACACTTCT'}='1c2';
$R{2}{'AAATGGCCGGAAAACAACACTTCT'}='1c3';
#$R{2}{'GAATGGCCGGAAGACAACACTTCT'}='1c4';
#$R{2}{'GAATGGCCCAAAGACAACGGCTCT'}='1c5';
$R{2}{'GAATGGCCCGCCGACAACGGCGCT'}='2c1';
$R{2}{'GAATGGCCGAAAGACAACGACAAG'}='2c2';
$R{2}{'ATATGGCCGGAAAACAACGACAAG'}='2c3';
$R{2}{'GAATGGCCGGAAGACAACACTTCT'}='1c4 2c4';
$R{2}{'GAATGGCCGGAAAACAACGACAAG'}='2c5';
#$R{2}{'ATATGGCCGAAAGACAACACTTCT'}='2c6';
$R{2}{'GAATGGCCGGAAAACAACACTTCT'}='3c1';
$R{2}{'GAATGGCCCAAAGACAACGGCTCT'}='1c5 3c2';
#$R{2}{'ATATGGCCGAAAGACAACACTTCT'}='3c3';
$R{2}{'AAATGGCCCGCCGACAACGGCGCT'}='6c1';
$R{2}{'AAATGGCCGGAAGACAACACTTCT'}='6c2';
$R{2}{'AAATGGCCCGCAAACAACGGCGCT'}='6c3';
$R{2}{'ACATGGCCGGAAAACAACACTTCT'}='7c1';
$R{2}{'ACATGGCCGAAAGACAACGGTGAT'}='1c1 uss';
$R{2}{'ATATGGCCGAAAGACAACACTTCT'}='ref'; # 2c6, 3c3
$R{2}{'ACATGGCCGAAAGACAACACTTCT'}='1c1 7c1 (T188A)';
$R{2}{'GAATGGCCGAAAGACAACACTTCT'}='var (A187G)';
$R{2}{'ATATGGCCGGAAAACAACACTTCT'}='var (195/199)';

# Region 3
$R{3}{'TCCCCCGCCGACAAAATCAA'}   = '1c1';
$R{3}{'TCCCCCGCCGAAATCAA'}      = '1c2';
$R{3}{'TCCCCCCCCCTCCGACATCAA'}  = '1c3';
$R{3}{'TCCCCCCCCACCGACATCAA'}   = '1c4';
$R{3}{'TCCGCTTCAAAAATCAT'}      = '1c5';
#$R{3}{'TCCGCTTCAACAATCAA'}      = '2c1';
$R{3}{'TCCCCCCCCTCCAACATCAA'}   = '2c2';
$R{3}{'TCTTCTTCATCAATCAA'}      = '2c3';
$R{3}{'TCCCCCCCACCGACATCAA'}    = '2c4';
#$R{3}{'TCCGCCTCCGACATCAA'}      = '2c5';
$R{3}{'TCCTCCGCCGCCGACATCAA'}   = '2c6';
$R{3}{'TCCTCCGACAAAATCAA'}      = '3c1';
#$R{3}{'TCCGCTTCAACAATCAA'}      = '3c2';
$R{3}{'AACCCCACCGACATCAA'}      = '3c3';
$R{3}{'TCCCCCGCCACCGACATCAA'}   = '6c1';
$R{3}{'TCCCCCCCTCCGACATCAA'}    = '6c2';
$R{3}{'TCCGCCTCCGACATCAA'}      = '2c5 6c3';
$R{3}{'TCCTCCGCCACCGACATCAA'}   = '7c1';
$R{3}{'TCCTCCGCCGAAATCAA'}      = 'uss';
$R{3}{'TCCGCTTCAACAATCAA'}      = 'ref'; #2c1, 3c2
$R{3}{'TCCGCTTCAACAATCAT'}      = '1c5 (A239T)';

# Region 4
#$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = '1c1';
#$R{4}{'AAAAGCGTTACGGTCGCAAAA'} = '1c2';
$R{4}{'AAAGAGGTTGAAGTTAAAAAC'} = '1c3';
#$R{4}{'AAAAGCGTTACGGTCGCAAAA'} = '1c4';
#$R{4}{'AAGGAAGTTAAAGTCGAAAAC'} = '1c5';
#$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = '2c1';
#$R{4}{'GAAAGCGTTACGGTCACAAAC'} = '2c2';
#$R{4}{'AAGGAAGTTAAAGTCGAAAAC'} = '2c3';
#$R{4}{'AAAAGCGTTACGGTCGCAAAA'} = '2c4';
$R{4}{'AAAAGCGTTACGGTCGCAAAA'} = '1c2 1c4 2c4 2c5'; 
$R{4}{'CAGAAAGTTGAAGTCAACAAC'} = '2c6';
#$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = '3c1';
$R{4}{'CAGAAAGTTGAAGTCACAAAC'} = '2c2 3c2 3c3 (268/273)';
$R{4}{'GAAAGCGTTACGGTCACAAAC'} = '2c2 3c3';
$R{4}{'AAGGAAGTTAAAGTCGAAAAC'} = '1c5 2c3 6c1';
#$R{4}{'CAAAGCGTTACGGTCGCAAAC'} = '6c2';
$R{4}{'GAAAGCGTTACGGTCGAAAAA'} = '6c3';
$R{4}{'CAAAGCGTTACGGTCGCAAAC'} = '6c2 7c1';
$R{4}{'AAAAGCGTTACGGTCGCAAAC'} = 'uss';
$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = 'ref'; #1c1, 2c1, 3c1
#CAGAAAGTTGAAGTCACAAAC = 2c2, 3c2, 3c3 (this is 268/273)
$R{4}{'CAGAAAGTTGAAGTCGAAAAC'} = '1c5 2c3 6c1 (269/273)';


# Region 5
#$R{5}{'AGGCGTCGTTACCGCCGA'} = '1c1';
#$R{5}{'AGGCGTCGTTACCGCCCA'} = '1c2';
#$R{5}{'CGGCGTCGTTACCGCCAC'} = '1c3';
#$R{5}{'AGGCGTCGTCACCGCCGA'} = '1c4';
#$R{5}{'CGGCGTCGTCACCGCCCA'} = '1c5';
#$R{5}{'AGGCGTCGTTACCGCCCA'} = '2c1';
#$R{5}{'CGGCGTCGTTACCGCCAC'} = '2c2';
$R{5}{'CGGCGTCGTCACCGCCAC'} = '2c3';
$R{5}{'AGGCGTCGTCACCGCCGA'} = '1c4 2c4';
$R{5}{'AGGCGTCGTTACCGCCGA'} = '1c1 2c5';
$R{5}{'CGGCGTCGTTACCGCCAC'} = '1c3 2c2 2c6';
#$R{5}{'AGGCGTCGTTACCGCCCA'} = '3c1';
$R{5}{'CGGCGTCGTTACCGCCCA'} = '3c2';
$R{5}{'CGGCGTCGTTACCGCCAA'} = '3c3';
$R{5}{'CGGCGTCGTCACCGCCCA'} = '1c5 6c1';
#$R{5}{'CGGCGTCGTTACCGCCGA'} = '6c2';
$R{5}{'AGGCGTCGTTACCGCCAA'} = '6c3';
#$R{5}{'CGGCGTCGTTACCGCCGA'} = '7c1';
$R{5}{'CGGCGTCGTTACCGCCGA'} = '6c2 7c1 uss';
$R{5}{'AGGCGTCGTTACCGCCCA'} = 'ref'; # 1c2, 2c1, 3c1
#HERE
# Region 6
$R{6}{'GAAATGAAACCAAGCGG'} = '1c1';
$R{6}{'CAAATGAATCCAAGCGG'} = '1c2';
#$R{6}{'ACAATGCTTTCAAGCGG'} = '1c3';
#$R{6}{'GAAATGGCTTCAACCGG'} = '1c4';
$R{6}{'CAAATGGCTTCAAGCAA'} = '1c5';
#$R{6}{'CAAATGGCTTCAACCGG'} = '2c1'; 
$R{6}{'ACAATGCTTTCAAGCGG'} = '1c3 2c2';
$R{6}{'ACAATGAATTCAAGCAA'} = '2c3';
#$R{6}{'GAAATGGCTTCAACCGG'} = '2c4';
$R{6}{'GAAATGGCTTCAACCGG'} = '1c4 2c4 2c5';
$R{6}{'ACAATGGCTTCAAGCAA'} = '2c6';
#$R{6}{'CAAATGGCTTCAACCGG'} = '3c1'; 
#$R{6}{'CAAATGGCTTCAACCGG'} = '3c2'; 
#$R{6}{'AAAATGCTTTCAAGCGG'} = '3c3';
#$R{6}{'CAAATGGCTTCAACCGG'} = '6c1'; 
#$R{6}{'GAAATGAAATCAGACGG'} = '6c2';
$R{6}{'AAAATGCTTTCAAGCGG'} = '3c3 6c3';
$R{6}{'GAAATGAAATCAGACGG'} = '6c2 7c1';
$R{6}{'GAAATGGCTTCAAGCGG'} = 'uss';
$R{6}{'CAAATGGCTTCAACCGG'} = 'ref'; # 2c1, 3c1, 3c2, 6c1

# Region 7
#$R{7}{'AGAAATCAAAGGCAAAAA'} = '1c1';
$R{7}{'TGAAATCAAAGACAAAAA'} = '1c2';
#$R{7}{'TGAAATCAAAGGCAAAAA'} = '1c3';
#$R{7}{'TGAAATCAAAGGCAAAAA'} = '1c4';
$R{7}{'AGAAATCAAAGACAAAAA'} = '1c5';
#$R{7}{'AGAAATCCAAGACAAAAA'} = '2c1';
#$R{7}{'TGAAATCAAAGGCAAAAA'} = '2c2';
$R{7}{'AGAAATCAAAGACAAAAG'} = '2c3';
#$R{7}{'TGAAATCAAAGGCAAAAA'} = '2c4';
#$R{7}{'AGAAATCCAAGGCAAAAA'} = '2c5';
#$R{7}{'AGAAATCCAAGGCAAAAG'} = '2c6';
#$R{7}{'AGAAATCCAAGGCAAAAA'} = '3c1';
#$R{7}{'AGAAATCAAAGGCAAAAA'} = '3c2';
#$R{7}{'AGAAATCCAAGGCAAAAG'} = '3c3';
$R{7}{'TGAAATCAAAGGCAAAAA'} = '1c3 1c4 2c2 2c4 6c1';
$R{7}{'AGAAATCAAAGGCAAAAA'} = '1c1 3c2 6c2 var322/326';
$R{7}{'AGAAATCCAAGGCAAAAA'} = '2c5 3c1 6c3 varA362G';
$R{7}{'AGAAATCCAAGGCAAAAG'} = '2c6 3c3 7c1';
$R{7}{'AGAAATCAAAGGCAAAAG'} = 'uss';
$R{7}{'AGAAATCCAAGACAAAAA'} = 'ref'; #2c1 
#AGAAATCCAAGGCAAAAA = var  (this is the A326G)
#AGAAATCAAAGGCAAAAA = var (this is 322/326)



# Region 8
$R{8}{'CCAAGCGTGAAGACGGTTC'} = '1c1';
#$R{8}{'CCAAGCGTGAAAACGGTTC'} = '1c2';
$R{8}{'CCAGGCGTGAAAACGGTTC'} = '1c3';
$R{8}{'CCAAGCGTCAAGACGGTC'}  = '1c4';
#$R{8}{'CCAAGCGTGAAAACGGTTC'} = '1c5';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '2c1';
$R{8}{'CCAAGCGTCAAGCCGGTTC'} = '2c2';
#$R{8}{'CCAAGCGTGAAAACGGTTC'} = '2c3';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '2c4';
$R{8}{'CCAAGCGTGAAAACGGTTC'} = '1c2 1c5 2c3 2c5 var355/358';
$R{8}{'CCAAGCGTGAAAACGGTC'}  = '2c6';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '3c1';
#$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2';
#$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '6c1';
$R{8}{'GCAGGCGTGAAAACGGTTC'} = '6c2';
$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3 6c3';
$R{8}{'CCAGGCGTGAAGCCGGTTC'} = '7c1';
$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2 uss varA350G';
$R{8}{'CCAAGCGTCAAGACGGTTC'} = 'ref'; # 2c1 2c4 3c1 6c1
#CCAGGCGTCAAGACGGTTC = var (this is the A350G)
#CCAAGCGTGAAAACGGTTC = var (this is 355/358)

# Region 9
$R{9}{'AAGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA'}                   = '1c1';
$R{9}{'ACGCGCAACGACGCCAAAGCCGACGCCAAAGACGACACCGTCACCGCCATCGAA'}                         = '1c2';
$R{9}{'ACGCGCGCCAAAGCCGACGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC'}                   = '1c3';
$R{9}{'AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAAAGACACCGACAAAATCAAC'}       = '1c4';
#$R{9}{'AAGCGCACCGAAGCCAACGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACACCGCCGGCACCAAAATCGAC'}    = '1c5';
#$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c1';
$R{9}{'CAGCGCGCCAAAGCCGACGACGCCGTCACCGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCGCCAACGTTGCCGCCGCCAACGACGACGACGTTACCGACGACAAAAACAACAACGGCATCGAC'}          = '2c3';
$R{9}{'AAGCGCGCCGACAACAACGGCAACATTACCGCCGACAACGGCAACGCCATCGAA'}                         = '2c4';
$R{9}{'AAGCGCACCGAAGCCAACGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACACCGCCGGCACCAAAATCGAC'}    = '1c5 2c5';
$R{9}{'ACGCGCGACGACAAAGCCAAAGACGACGTCAAAGCCGACGGCACCGCCGGCACCAAAATCGAC'}                = '2c6';
$R{9}{'ACGCGCAACGACGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA'}                      = '3c1';
$R{9}{'AAGCGCGACGACGCCGCCGCCAAAGACGACACCGTCACCGCCGACGCCACCGGCAACGACGGCAAAATCGAC'}       = '3c2';
$R{9}{'AAGCGCACCGAAGCCAACGCCGACGCCGCCGGCAAAGACACCACCAACGGCATCAAC'}                      = '3c3';
$R{9}{'AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAACAACGGCATCAAC'}             = '6c1';
$R{9}{'AAGCGCGACGCCAACAACGCCAACAACGACGCCGTCACCGACGACACCACCGGCAACGGCAACGAAAAAATCGAA'}    = '6c2';
$R{9}{'AAGCGCAACGACGCCGCCAACGACGACGTTACCGACGACGCCGGCACCGACAACGGCGGCAAAGGCAAAATCGAC'}    = '6c3';
$R{9}{'ACGCGCGCCAAAGCCAAAGACGCCGACGACGTTACCGACGACGCCGGCACCCACAACGGCGGCAAAGGCAAAATCGAC'} = '7c1';
$R{9}{'ACGCGCAACGACGCCGCCGACAACGACGACGTCGCCAAAGACGACGCCGCCGGCAACGCCATCGAA'}             = 'uss';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'ref'; #2c1
$R{9}{'CAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'var';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAA'}                         = 'var';

 
#TODO
#Ccgccaac = var (this is 431 – 438)
#R9.5 (445 – 452)
#GACACCAGG = 1c5 (this is A452G)

# Region 10
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '1c1';
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '1c2';
$R{10}{'CGATGAATCATCTGCCGTT'}    = '1c3';
$R{10}{'CGATAAATCATCTGCCGTT'}    = '1c4';
#$R{10}{'CGATGAATCATCGTTGCCGG'}   = '1c5';
$R{10}{'CGATAAACATGATGCCAAATG'}  = '2c1';
$R{10}{'TGATACGTCATCTGCCAAA'}    = '2c2';
$R{10}{'CGATAAATCATCTGCCACCTA'}  = '2c3 var';
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '2c4';
#$R{10}{'CGATGAATCATCGTTGCCGG'}   = '2c5';
$R{10}{'CGATGAATCATCGTTGCCGG'}   = '1c5 2c5 2c6';
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '3c1';
#$R{10}{'CGATAAATCAACTGCCGTT'}    = '3c2';
$R{10}{'CGACCCGTTCTCTGCTAGC'}    = '3c3';
$R{10}{'CGATAAACATGATGCCAAA'}    = '6c1';
$R{10}{'CGATGAATCATCTGCCGTTTA'}  = '6c2';
$R{10}{'CGATAAATCAACTGCCGTT'}    = '3c2 6c3';
$R{10}{'CGATAAATCAACTGCCAAA'}    = '7c1';
$R{10}{'CGATGAACCAACTGCCACCTA'}  = 'uss';
$R{10}{'CGATGAATCATCTGCCACCTA'}  = 'ref'; #1c1 1c2 2c4 3c1
$R{10}{'TGATGAATCATCTGCCACCTA'}  = '2c2';
$R{10}{'CGATGAATCATCTGCCAAATA'}  = 'var';
#CGATAAATCATCTGCCACCTA = var


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
    #TODO continue testing the latter part of the sequence
    next;
  }

  for ($r=1; $r<=$nRegions; $r++){ #start later for debugging
    while (my ($key_seq,$silent_copy) = each $R{$r}){
    
      if ($debug>=1) {
	      print "\$r=$r, $silent_copy\n"
      }
      my $substr = substr $seq_obj->seq, $region_min[$r]-$offset-$fudge_factor, $region_length[$r]+2*$fudge_factor;
      my $search = $key_seq;
      if ($debug>=1) {
	      print "\$substr=".$substr."\n";
	      print "\$search=".$search."\n";
      }
      if ( $substr =~ /$search/ ){
	      print $silent_copy," ";
        #if ($silent_copy eq "ref"){
        #   last;  # if it's ref, don't bother checking the others
        #}
      } 
    } 
        print ",";
  }
  print "\n";


}
