#!/usr/bin/perl -w

# find key sequences in an input file.
# a keyword for incomplete code is TODO

# To run: ./SilentDNA.pl inputfile.fna [ first_read ] [ last_read ] > outputfile.csv

# sequences are labeled according to specific silent copy names.
# if a sequence is tagged but has no specific silent copy, the keyword is "var"

# DONE   tolerate short reads
# DELETE implement short regions in region 9 
# TODO   implement a "verdict" column
# TODO   try to index region 10 from the right, since R9 length varies
# TODO   count the number of empty cells on each row
# TODO   columns for counts of "var (any incl named)", "blank", and "ref"

use Bio::SeqIO;

#modify the inputfile filename as needed
#$inputfile = "SampleFile.txt";
#$inputfile = "MID24_FA1090.fna";
#$inputfile = "MID4_FA1090.fna";

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

$debug=0; #set to 1 for some debugging output

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
    290-273+1,
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
$ref_length[1]=length('CGTCACCGAGTATTACCTGAATA');

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
$ref_length[2]=length('ATATGGCCGAAAGACAACACTTCT');

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
$ref_length[4]=length('CAGAAAGTTGAAGTCGCAAAA');


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
$ref_length[5]=length('AGGCGTCGTTACCGCCCA');

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
$ref_length[6]=length('CAAATGGCTTCAACCGG');

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
$ref_length[7]=length('AGAAATCCAAGACAAAAA');



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
$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3 6c3';
$R{8}{'CCAGGCGTGAAGCCGGTTC'} = '7c1';
#$R{8}{'GCGTGAAGCCGGTTC'} = 'v355/359';
$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2 uss vA350G';
$R{8}{'CCAAGCGTCAAGACGGTTC'} = 'ref'; # 2c1 2c4 3c1 6c1
#CCAGGCGTCAAGACGGTTC = var (this is the A350G)
#CCAAGCGTGAAAACGGTTC = var (this is 355/358)
$ref_length[8]=length('CCAAGCGTCAAGACGGTTC');

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
$R{9}{'ACGCGCAACGACGCCGCCGACAACGACGACGTCGCCAAAGACGACGCCGCCGGCAACGCCATCGAA'}             = 'uss';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'ref'; #2c1
$R{9}{'CAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'var';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAA'}                         = 'var';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCGCCAACGCCATCGAC'}                          = 'var';  # added Oct 15, 2013
$ref_length[9]=length('ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC');

 
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
$R{10}{'CGATAAATCATCTGCCACCTA'}  = '2c3'; #also a var
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
$R{10}{'CGATGAATCAACTGCCAAA'}    = '7c1'; # added Oct 15, 2013
$R{10}{'CGATGAACCAACTGCCACCTA'}  = 'uss';
$R{10}{'CGATGAATCATCTGCCACCTA'}  = 'ref'; #1c1 1c2 2c4 3c1
$R{10}{'CGATGAATCATCTGCCAAATA'}  = 'var491';
$R{10}{'CGATGAATCATCTGTCACCTA'}  = 'false';
$ref_length[10]=length('CGATGAATCATCTGCCACCTA');

$bound_09 = "GGTAAAATGGTTCTGCGGACAGCCGGTT"; #boundary between regions 9 and 10


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

sub printVerdict{
    my @result = @_;
    my $resultString = join(",",@result);
    #print $resultString."\n";
    my $nRef =()= $resultString =~ /ref/gi;
    if ($nRef==10){
        #if all are ref, print ref
        print "ref,";
    } elsif ($resultString =~ /var/){
        #if it's var, just print var
        print "var,";
    } elsif (getNblank(@result) != 0){
        print ","; #blank
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
            print $reversebitmap{$totalBitString}.",";
        } else {
            print "var,";
        }
    }

}



#header
print "counter,ID,region1,region2,region3,region4,region5,region6,region7,region8,region9,region10,nVar,nEmpty,verdict,basepairs\n";
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
    #check sequence is long enough to run the tests
    if (length($seq_obj->seq) < $region_min[$nRegions]+$region_length[$nRegions]){ #+2*$fudge_factor){
        #print "error: sequence too short\n";
        $this_offset=271;
        print "short,short,short,short,short,";
        $rStart=6;
        for ($r=1; $r<6; $r++){
            $result[$r]="short";
        }
        #TODO continue testing the latter part of the sequence
        #next; #go to next sequence
    }


    for ($r=$rStart; $r<=$nRegions; $r++){ 
        my $substr;
        if ( ($r==10) && ($debug>=1) )#index region 10 from the end
        {
            #find the match position
            $seq_obj->seq =~ /$bound_09/;
            print "region 9 to 10 boundary found between @- and @+ \n";
        } 
        $substr = substr $seq_obj->seq, $region_min[$r]-$this_offset-$fudge_factor+$length_shift, $region_length[$r]+2*$fudge_factor;


        while (my ($key_seq,$silent_copy) = each $R{$r}){

            if ($debug>=1) {
                print "\$r=$r, $silent_copy\n"
            }
            # substr source start length
 
            my $search = $key_seq;
            if ($debug>=1) {
                print "\$substr=".$substr."\n";
                print "\$search=".$search."\n";
            }
            if ( $substr =~ /$search/ ){
                print $silent_copy," ";
                if ( ($r==10) && ($debug>=1) ){print "region 10 found between @- and @+ \n";}
                $result[$r]=$silent_copy;
                $length_shift += length($search) - $ref_length[$r];
                if ($debug>=1) {
                    print @-,"-",@+," ";
                    print $length_shift;
                }
                #if ($silent_copy eq "ref"){
                #   last;  # if it's ref, don't bother checking the others
                #}
            }        
        } 
        if ($result[$r] eq " "){
          print $substr; # if no match found, print the search sequence to figure out why
        }
        print ",";
    }

    &printNvar(@result); # print the number of variants
    &printNblank(@result); # print the number of blanks
    &printVerdict(@result); # print the verdict

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
