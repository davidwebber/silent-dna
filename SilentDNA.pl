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

$offset=91;  # offset between reference sequence and reads

$fudge_factor = 3;  #the number of insertions or deletions allowed before the sequence of interest


##########################################################
# you shouldn't have to modify anything below this line
##########################################################


#$nSilent=20; #number of silent copies, including the reference

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
    288-274+1,
    305-289+1,
    332-315+1,
    365-347+1,
    447-394+25,
    494-474+1);#}}}

#@index=("ref", #{{{
#        "1c1",
#        "1c2",
#        "1c3",
#        "1c4",
#        "1c5",
#        "2c1",
#        "2c2",
#        "2c3",
#        "2c4",
#        "2c5",
#        "2c6",
#        "3c1",
#        "3c2",
#        "3c3",
#        "6c1",
#        "6c2",
#        "6c3",
#        "7c1",
#        "uss"); #}}}

# Region 1 
#$R{1}{'CGTTGCCGGGTATTGCCCGAATC'}='1c1';
#$R{1}{'CGTTACCGGGTATTGCCCGAATC'}='1c2';
$R{1}{'CGTTGCCGGGTATTGCCCGAATC'}='1c1 1c3';
$R{1}{'AGTTGCCGGGTATTGCCTGAATC'}='1c4';
#$R{1}{'CGTTACCGAGTATTACCTGAATC'}='1c5';
$R{1}{'CGTCACCGAATATTACCCGAATA'}='2c1';
$R{1}{'CGTTACCGAGTATTGCCCGAATC'}='2c2';
$R{1}{'CGTTACCGAGTATTACCTGAATC'}='1c5 2c3 v162';
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
$R{1}{'CGTTGCCGAGTATTACCTGAATC'}='v163';
$R{1}{'CGTCACCGAGTATTACCTGAATA'}='2c1 6c1 6c3 (C181A)';
#$R{1}{'CGTCACCGAGTATTACCTGAAAA'} = 'indel'; #same as 3c2 hybrid
#$R{1}{'CGTCACCGAGTATTACCCGAATC'} = 'var (T176C)'; #same as 6c2 hybrid
$R{1}{'CGTCACCGAGTATTACCCGAATC'} = '6c2 hybrid';
$R{1}{'CGTCACCGGGTATTACCTGAATA'} = '6c1 6c3 2c1 hybrid';
$R{1}{'CGTCACCGAGTATTACCTGAAAA'} = '3c2 hybrid';
$ref_length[1]=length('CGTCACCGAGTATTACCTGAATA');
$R{1}{'ref'}='CGTCACCGAGTATTACCTGAATC'; # uss, 7c1, 3c3


# Region 2
#$R{2}{'ACATGGCCGAAAGACAACGGTGAT'}='1c1';
$R{2}{'ATATGGCCGAAAGACAACGGTGAT'}='1c1'; # added Oct 15, 2013
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
$R{2}{'GAATGGCCGAAAGAACAACACTTCT'}='var (SF37)';  # added Oct 15, 2013
$R{2}{'ATATGGCCGAAAGACAACACGCT'} ='indel';
$R{2}{'ATATGGACCGAAAGACAACACTTCT'}='indel';
$R{2}{'ATATGGTCGAAAGACAACACTTCT'} = 'false';
$R{2}{'ATATGGCCCGCCGACAACGGCGCT'} = '2c1 6c1 hybrid';
$ref_length[2]=length('ATATGGCCGAAAGACAACACTTCT');
$R{2}{'ref'}='ATATGGCCGAAAGACAACACTTCT'; # 2c6, 3c3

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
$R{3}{'TCCGCGACAACAATCAA'}      = 'false';
$R{3}{'TCCGACTTCAACAATCAA'}     = 'indel';
$ref_length[3]=length('TCCGCTTCAACAATCAA');
$R{3}{'ref'}='TCCGCTTCAACAATCAA'; #2c1, 3c2

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
$R{4}{'CAGAAAGTTGAAGTCACAAAC'} = '2c2 3c2 3c3 v(268/273)';
$R{4}{'GAAAGCGTTACGGTCACAAAC'} = '2c2 3c3';
$R{4}{'AAGGAAGTTAAAGTCGAAAAC'} = '1c5 2c3 6c1';
#$R{4}{'CAAAGCGTTACGGTCGCAAAC'} = '6c2';
$R{4}{'GAAAGCGTTACGGTCGAAAAA'} = '6c3';
$R{4}{'CAAAGCGTTACGGTCGCAAAC'} = '6c2 7c1';
$R{4}{'AAAAGCGTTACGGTCGCAAAC'} = 'uss';
$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = 'ref'; #1c1, 2c1, 3c1
#CAGAAAGTTGAAGTCACAAAC = 2c2, 3c2, 3c3 (this is 268/273)
$R{4}{'CAGAAAGTTGAAGTCGAAAAC'} = '1c5 2c3 6c1 v(269/273)';
$R{4}{'CAGAAAGTTGAAGTCGCAAAG'} = 'indel';
$R{4}{'AAGGAAGTTAAAGTCGCAAAA'} = '1c5 2c3 6c1 hybrid';
$R{4}{'CAGAAAGTTGAAGTCGCAAAC'} = '6c2 7c1 uss hybrid';
$R{4}{'CAGAAAGTTGAAGTCACAAAA'} = '3c2 3c3 2c2 hybrid';
$ref_length[4]=length('CAGAAAGTTGAAGTCGCAAAA');
$R{4}{'ref'}='CAGAAAGTTGAAGTCGCAAAA'; #1c1, 2c1, 3c1

## Region 5
##$R{5}{'AGGCGTCGTTACCGCCGA'} = '1c1';
##$R{5}{'AGGCGTCGTTACCGCCCA'} = '1c2';
##$R{5}{'CGGCGTCGTTACCGCCAC'} = '1c3';
##$R{5}{'AGGCGTCGTCACCGCCGA'} = '1c4';
##$R{5}{'CGGCGTCGTCACCGCCCA'} = '1c5';
##$R{5}{'AGGCGTCGTTACCGCCCA'} = '2c1';
##$R{5}{'CGGCGTCGTTACCGCCAC'} = '2c2';
#$R{5}{'CGGCGTCGTCACCGCCAC'} = '2c3';
#$R{5}{'AGGCGTCGTCACCGCCGA'} = '1c4 2c4';
#$R{5}{'AGGCGTCGTTACCGCCGA'} = '1c1 2c5';
#$R{5}{'CGGCGTCGTTACCGCCAC'} = '1c3 2c2 2c6';
##$R{5}{'AGGCGTCGTTACCGCCCA'} = '3c1';
#$R{5}{'CGGCGTCGTTACCGCCCA'} = '3c2';
#$R{5}{'CGGCGTCGTTACCGCCAA'} = '3c3';
#$R{5}{'CGGCGTCGTCACCGCCCA'} = '1c5 6c1';
##$R{5}{'CGGCGTCGTTACCGCCGA'} = '6c2';
#$R{5}{'AGGCGTCGTTACCGCCAA'} = '6c3';
##$R{5}{'CGGCGTCGTTACCGCCGA'} = '7c1';
#$R{5}{'CGGCGTCGTTACCGCCGA'} = '6c2 7c1 uss';
#$R{5}{'AGGCGTCGTTACCGCCCA'} = 'ref'; # 1c2, 2c1, 3c1
#$ref_length[5]=length('AGGCGTCGTTACCGCCCA');

# Region 5
$R{5}{'GGCGTCGTTACCGCC'} = 'ref'; # 1c2, 2c1, 3c1
$R{5}{'GGCGTCGTCACCGCC'} = 'v282 1c4 1c5 2c3 2c4 6c1';
$ref_length[5]=length('GGCGTCGTTACCGCC');
$R{5}{'ref'}='GGCGTCGTTACCGCC'; # 1c2, 2c1, 3c1

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
$R{6}{'CAAATGACTTCAACCGG'} = 'false';
$R{6}{'CAAATGGCTTCAAACCGG'} = 'indel';
$ref_length[6]=length('CAAATGGCTTCAACCGG');
$R{6}{'ref'}='CAAATGGCTTCAACCGG'; # 2c1, 3c1, 3c2, 6c1

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
$R{7}{'AGAAATCAAAGGCAAAAA'} = '1c1 3c2 6c2 v322/326';
$R{7}{'AGAAATCCAAGGCAAAAA'} = '2c5 3c1 6c3 vA362G';
$R{7}{'AGAAATCCAAGGCAAAAG'} = '2c6 3c3 7c1';
$R{7}{'AGAAATCAAAGGCAAAAG'} = 'uss';
$R{7}{'AGAAATCCAAGACAAAAA'} = 'ref'; #2c1 
#AGAAATCCAAGGCAAAAA = var  (this is the A326G)
#AGAAATCAAAGGCAAAAA = var (this is 322/326)
$R{7}{'AGAAATACCAAGACAAAAA'} = 'indel';
$R{7}{'AGAAATCCAAGCAAAAA'} = 'indel';
$R{7}{'AGAAATCCAAAGACAAAAA'} = 'indel';
$R{7}{'AGAAAATCCAAGACAAAAA'} = 'indel';
$R{7}{'TGAAATCCAAGACAAAAA'} = 'var (A315T)';
$R{7}{'AGAAATCCAAAGGCAAAAGA'} = '2c6 3c3 7c1 hybrid';
$ref_length[7]=length('AGAAATCCAAGACAAAAA');
$R{7}{'ref'}='AGAAATCCAAGACAAAAA'; #2c1 



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
$R{8}{'CCAAGCGTGAAAACGGTTC'} = '1c2 1c5 2c3 2c5 v355/358';
$R{8}{'CCAAGCGTGAAAACGGTC'}  = '2c6';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '3c1';
#$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2';
#$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3';
#$R{8}{'CCAAGCGTCAAGACGGTTC'} = '6c1';
$R{8}{'GCAGGCGTGAAAACGGTTC'} = '6c2';
#$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3 6c3';
$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3 6c3 hybrid 7c1';
$R{8}{'CCAGGCGTGAAGCCGGTTC'} = '7c1';
#$R{8}{'GCGTGAAGCCGGTTC'} = 'v355/359';
$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2 uss vA350G';
$R{8}{'CCAAGCGTCAAGACGGTTC'} = 'ref'; # 2c1 2c4 3c1 6c1
#CCAGGCGTCAAGACGGTTC = var (this is the A350G)
#CCAAGCGTGAAAACGGTTC = var (this is 355/358)
$R{8}{'CCAAAGCGTCAAGACGGTTC'} = 'indel';
$ref_length[8]=length('CCAAGCGTCAAGACGGTTC');
$R{8}{'ref'}='CCAAGCGTCAAGACGGTTC'; # 2c1 2c4 3c1 6c1

# Region 9
$R{9}{'AAGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA'}                   = '1c1';
$R{9}{'ACGCGCAACGACGCCAAAGCCGACGCCAAAGACGACACCGTCACCGCCATCGAA'}                         = '1c2';
$R{9}{'ACGCGCGCCAAAGCCGACGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC'}                   = '1c3';
$R{9}{'AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAAAGACACCGACAAAATCAAC'}       = '1c4';
#$R{9}{'AAGCGCACCGAAGCCAACGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACACCGCCGGCACCAAAATCGAC'}    = '1c5';
#$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c1';
$R{9}{'CAGCGCGCCAAAGCCGACGACGCCGTCACCGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCGCCAACGTTGCCGCCGCCAACGACGACGACGTTACCGACGACAAAAACAACAACGGCATCGAC'}          = '2c3';
$R{9}{'AGCGCGCCAACGTTGCCGCCGCCAACGACGACGACGTTACCGACGCCAACAACGCCATCGAC'}                 = '2c3';  # added Oct 15, 2013
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
$R{9}{'ACGCGCGCCAAAGCCAAAGACGCCGACGACGTTACCGACGACGCCGGCACCGACAACGGCGGCAAAGGCAAAATCGAC'} = '7c1';
$R{9}{'ACGCGCAACGACGCCGCCGACAACGACGACGTCGCCAAAGACGACGCCGCCGGCAACGCCATCGAA'}             = 'uss';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'ref'; #2c1
$R{9}{'CAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'var';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAA'}                         = 'var';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCGCCAACGCCATCGAC'}                          = 'var';  # added Oct 15, 2013
$ref_length[9]=length('ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC');
$R{9}{'ref'}='ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'; #2c1

 
#TODO
#Ccgccaac = var (this is 431 – 438) # Ella will send me a longer sequence
#R9.5 (445 – 452)
#GACACCAGG = 1c5 (this is A452G) # run this as a special case

# Region 10
# If the regions are shorter than the reference, then they just need a "TA" at the end
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '1c1';
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '1c2';
$R{10}{'CGATGAATCATCTGCCGTT'}    = '1c3';
$R{10}{'CGATAAATCATCTGCCGTT'}    = '1c4';
#$R{10}{'CGATGAATCATCGTTGCCGG'}   = '1c5';
#$R{10}{'CGATAAACATGATGCCAAATG'}  = '2c1';
$R{10}{'TGATACGTCATCTGCCAAA'}    = '2c2';
$R{10}{'TGATACGTCATCTGCCACCTA'}  = '2c2'; # added Oct 15, 2013
$R{10}{'CGATAAATCATCTGCCACCTA'}  = '2c3 hybrid 7c1'; #also a var
#CGATAAATCATCTGCCACCTA = var
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '2c4';
#$R{10}{'CGATGAATCATCGTTGCCGG'}   = '2c5';
$R{10}{'CGATGAATCATCGTTGCCGG'}   = '1c5 2c5 2c6';
#$R{10}{'CGATGAATCATCTGCCACCTA'}  = '3c1';
#$R{10}{'CGATAAATCAACTGCCGTT'}    = '3c2';
$R{10}{'CGACCCGTTCTCTGCTAGC'}    = '3c3';
$R{10}{'CGATAAACATGATGCCAAA'}    = '2c1 6c1'; #2c1 has "TG" appended
$R{10}{'CGATGAATCATCTGCCGTTTA'}  = '6c2';
$R{10}{'CGATAAATCAACTGCCGTT'}    = '3c2 6c3';
$R{10}{'CGATAAATCAACTGCCAAA'}    = '7c1';
$R{10}{'CGATAAAATCAACTGCCAAA'}   = '7c1';
$R{10}{'CGATGAATCAACTGCCAAA'}    = '7c1'; # added Oct 15, 2013
$R{10}{'CGATGAACCAACTGCCACCTA'}  = 'uss';
$R{10}{'CGATGAATCATCTGCCACCTA'}  = 'ref'; #1c1 1c2 2c4 3c1
$R{10}{'CGATGAATCATCTGCCAAATA'}  = 'var491';
$R{10}{'CGATGAATCATCTGTCACCTA'}  = 'false';
$R{10}{'CGATGAATCATCTGCCACCCTA'} = 'indel';
$R{10}{'CGATGAAATCATCTGCCACCTA'} = 'indel';
$R{10}{'CGATGAATCATCTGCCCACCTA'} = 'indel';
$R{10}{'CGATAAATCAAACTGCCAAATA'} = '7c1';
$ref_length[10]=length('CGATGAATCATCTGCCACCTA');
$R{10}{'ref'}='CGATGAATCATCTGCCACCTA'; #1c1 1c2 2c4 3c1

$bound[2]  = "ACGGC";           #boundary between regions 1 and 2
$bound[3]  = "GCCGGCGTGGCA";    #boundary between regions 2 and 3
$bound[4]  = "AGGCAAATATGTT";   #boundary between regions 3 and 4
$bound[7]  = "CGTAAACAA";       #boundary between regions 6 and 7
$bound[8]  = "CTCTCCCTGTGGG";   #boundary between regions 7 and 8
$bound[9]  = "GGTAAAATGGTTCTGCGGACAGCCGGTT";  #boundary between regions 8 and 9
$bound[10] = "ACCAAGCACCTGCCGTCAACCTGCCG";    #boundary between regions 9 and 10
#$bound_10_1c5 = "ACCAGGCACCTGCCGTCAACCTGCCG"; #boundary between regions 9 and 10 specifically for 1c5

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

#$bitmap{'var'}   = 0b0000000000000000000;
$bitmap{'1c1'}   = 0b0000000000000000001;
$bitmap{'1c2'}   = 0b0000000000000000010;
$bitmap{'1c3'}   = 0b0000000000000000100;
$bitmap{'1c4'}   = 0b0000000000000001000;
$bitmap{'1c5'}   = 0b0000000000000010000;
$bitmap{'2c1'}   = 0b0000000000000100000;
$bitmap{'2c2'}   = 0b0000000000001000000;
$bitmap{'2c3'}   = 0b0000000000010000000;
$bitmap{'2c4'}   = 0b0000000000100000000;
$bitmap{'2c5'}   = 0b0000000001000000000;
$bitmap{'2c6'}   = 0b0000000010000000000;
$bitmap{'3c1'}   = 0b0000000100000000000;
$bitmap{'3c2'}   = 0b0000001000000000000;
$bitmap{'3c3'}   = 0b0000010000000000000;
$bitmap{'6c1'}   = 0b0000100000000000000;
$bitmap{'6c2'}   = 0b0001000000000000000;
$bitmap{'6c3'}   = 0b0010000000000000000;
$bitmap{'7c1'}   = 0b0100000000000000000;
$bitmap{'uss'}   = 0b1000000000000000000;
$bitmap{'ref'}   = 0b1111111111111111111;
#$bitmap{' '}     = 0b1111111111111111111;
#$bitmap{'short'} = 0b1111111111111111111;

while (($key, $value) = each %bitmap){
    $reversebitmap{$value} = $key;
}
$reversebitmap{0b0000000000000000000} = 'double_cross';
$reversebitmap{0b0000100000000100000} = '2c1_6c1';

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
                $length_shift = "@+" - ($region_min[$r]-$this_offset-1); #assignment, not modification;
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
