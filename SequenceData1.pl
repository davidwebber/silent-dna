#This sequence data for FA1090

# offset between reference sequence and reads
$offset=91; 

$nRegions=10;

@region_min=(0, 
    159,
    187,
    223,
    253,
    273,
    289,
    315,
    347,
    394,
    474);

@region_length=(0, 
    181-159+1, #make the computer do the subtraction
    210-187+1,
    239-223+4,
    273-253+1,    
    288-274+1,
    305-289+1,
    332-315+1,
    365-347+1,
    447-394+25,
    494-474+1);

$bound[2]  = "ACGGC";           #boundary between regions 1 and 2, original
$bound[3]  = "GCCGGCGTGGCA";    #boundary between regions 2 and 3
$bound[4]  = "AGGCAAATATGTT";   #boundary between regions 3 and 4
$bound[7]  = "CGTAAACAA";       #boundary between regions 6 and 7
$bound[8]  = "CTCTCCCTGTGGG";   #boundary between regions 7 and 8
$bound[9]  = "GGTAAAATGGTTCTGCGGACAGCCGGTT";  #boundary between regions 8 and 9
$bound[10] = "ACCA[AG]GCACCTGCCGTCAACCTGCCG";    #boundary between regions 9 and 10 for either 1c5 or ref

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

while (($key, $value) = each %bitmap){
    $reversebitmap{$value} = $key;
}
$reversebitmap{0b0000000000000000000} = 'double_cross';
$reversebitmap{0b0000100000000100000} = '2c1_6c1';

# Region 1 
$R{1}{'CGTTGCCGGGTATTGCCCGAATC'}='1c1 1c3';
$R{1}{'AGTTGCCGGGTATTGCCTGAATC'}='1c4';
$R{1}{'CGTCACCGAATATTACCCGAATA'}='2c1';
$R{1}{'CGTTACCGAGTATTGCCCGAATC'}='2c2';
$R{1}{'CGTTACCGAGTATTACCTGAATC'}='1c5 2c3 v162';
$R{1}{'GGTTGCCGGGTATTGCCTGAATC'}='2c4';
$R{1}{'CGTTACCGGGTATTGCCCGAATC'}='1c2 2c6';
$R{1}{'CGTTGCCGGGTATTACCTGAATC'}='2c5 3c1';
$R{1}{'CCTCACCGAATATTACCTGAAA' }='3c2';
$R{1}{'CGTTACCGGGTATTACCTGAATA'}='6c1';
$R{1}{'CGTTGCCGGGTATTACCCGAATC'}='6c2';
$R{1}{'CGTCACCGAGTATTACCCGAATA'}='2c1 6c3';
$R{1}{'CGTCACCGAGTATTACCTGAATC'}='ref 3c3 7c1'; # uss, 7c1, 3c3
$R{1}{'CGTTGCCGAGTATTACCTGAATC'}='v163';
$R{1}{'CGTCACCGAGTATTACCTGAATA'}='2c1 6c1 6c3 (C181A)';
$R{1}{'CGTCACCGAGTATTACCCGAATC'} = '6c2 hybrid';
$R{1}{'CGTCACCGGGTATTACCTGAATA'} = '6c1 6c3 2c1 hybrid';
$R{1}{'CGTCACCGAGTATTACCTGAAAA'} = '3c2 hybrid';
$ref_length[1]=length('CGTCACCGAGTATTACCTGAATA');
$R{1}{'ref'}='CGTCACCGAGTATTACCTGAATC'; # uss, 7c1, 3c3
$R{1}{'CGTTACCGGGTATTACCTGAATC'} = 'var';

# Region 2
$R{2}{'ATATGGCCGAAAGACAACGGTGAT'}='1c1'; # added Oct 15, 2013
$R{2}{'AAATGGCCGAAAGACAACACTTCT'}='1c2';
$R{2}{'AAATGGCCGGAAAACAACACTTCT'}='1c3';
$R{2}{'GAATGGCCCGCCGACAACGGCGCT'}='2c1';
$R{2}{'GAATGGCCGAAAGACAACGACAAG'}='2c2';
$R{2}{'ATATGGCCGGAAAACAACGACAAG'}='2c3';
$R{2}{'GAATGGCCGGAAGACAACACTTCT'}='1c4 2c4';
$R{2}{'GAATGGCCGGAAAACAACGACAAG'}='2c5';
$R{2}{'GAATGGCCGGAAAACAACACTTCT'}='3c1';
$R{2}{'GAATGGCCCAAAGACAACGGCTCT'}='1c5 3c2';
$R{2}{'AAATGGCCCGCCGACAACGGCGCT'}='6c1';
$R{2}{'AAATGGCCGGAAGACAACACTTCT'}='6c2';
$R{2}{'AAATGGCCCGCAAACAACGGCGCT'}='6c3';
$R{2}{'ACATGGCCGGAAAACAACACTTCT'}='7c1';
$R{2}{'ACATGGCCGAAAGACAACGGTGAT'}='1c1 uss';
$R{2}{'ATATGGCCGAAAGACAACACTTCT'}='ref 2c6 3c3'; # 2c6, 3c3
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
$R{2}{'TATGGCCCAAAGACAACGGCTCTG'}= '1c5 hybrid';
$R{2}{'TATGGCCGAAAGACAACGGCTCTG'}= 'var';
$R{2}{'GCCCGGATCAACCCGGGCGGCTTG'}= 'GitHub';

# Region 3
$R{3}{'TCCCCCGCCGACAAAATCAA'}   = '1c1';
$R{3}{'TCCCCCGCCGAAATCAA'}      = '1c2';
$R{3}{'TCCCCCCCCCTCCGACATCAA'}  = '1c3';
$R{3}{'TCCCCCCCCACCGACATCAA'}   = '1c4';
$R{3}{'TCCGCTTCAAAAATCAT'}      = '1c5';
$R{3}{'TCCGCTTCAAAATCAT'}       = '1c5 indel';
$R{3}{'TCCGCTTCAACAATCAT'}      = '1c5 (A239T)';
$R{3}{'TCCCCCCCCTCCAACATCAA'}   = '2c2';
$R{3}{'TCTTCTTCATCAATCAA'}      = '2c3';
$R{3}{'TCCCCCCCACCGACATCAA'}    = '2c4';
$R{3}{'TCCTCCGCCGCCGACATCAA'}   = '2c6';
$R{3}{'TCCTCCGACAAAATCAA'}      = '3c1';
$R{3}{'AACCCCACCGACATCAA'}      = '3c3';
$R{3}{'TCCCCCGCCACCGACATCAA'}   = '6c1';
$R{3}{'TCCCCCCCTCCGACATCAA'}    = '6c2';
$R{3}{'TCCGCCTCCGACATCAA'}      = '2c5 6c3';
$R{3}{'TCCTCCGCCACCGACATCAA'}   = '7c1';
$R{3}{'TCCTCCGCCGAAATCAA'}      = 'uss';
$R{3}{'TCCGCTTCAACAATCAA'}      = 'ref 2c1 3c2'; #2c1, 3c2
$R{3}{'TCCGCGACAACAATCAA'}      = 'false';
$R{3}{'TCCGACTTCAACAATCAA'}     = 'indel';
$ref_length[3]=length('TCCGCTTCAACAATCAA');
$R{3}{'ref'}='TCCGCTTCAACAATCAA'; #2c1, 3c2


# Region 4
$R{4}{'AAAGAGGTTGAAGTTAAAAAC'} = '1c3';
$R{4}{'AAAAGCGTTACGGTCGCAAAA'} = '1c2 1c4 2c4 2c5'; 
$R{4}{'CAGAAAGTTGAAGTCAACAAC'} = '2c6';
$R{4}{'CAGAAAGTTGAAGTCACAAAC'} = '2c2 3c2 3c3 v(268/273)';
$R{4}{'GAAAGCGTTACGGTCACAAAC'} = '2c2 3c3';
$R{4}{'AAGGAAGTTAAAGTCGAAAAC'} = '1c5 2c3 6c1';
$R{4}{'GAAAGCGTTACGGTCGAAAAA'} = '6c3';
$R{4}{'CAAAGCGTTACGGTCGCAAAC'} = '6c2 7c1';
$R{4}{'AAAAGCGTTACGGTCGCAAAC'} = 'uss';
$R{4}{'CAGAAAGTTGAAGTCGCAAAA'} = 'ref 1c1 2c1 3c1'; #1c1, 2c1, 3c1
$R{4}{'CAGAAAGTTGAAGTCGAAAAC'} = '1c5 2c3 6c1 v(269/273)';
$R{4}{'CAGAAAGTTGAAGTCGCAAAG'} = 'indel';
$R{4}{'AAGGAAGTTAAAGTCGCAAAA'} = '1c5 2c3 6c1 hybrid';
$R{4}{'CAGAAAGTTGAAGTCGCAAAC'} = '6c2 7c1 uss hybrid';
$R{4}{'CAGAAAGTTGAAGTCACAAAA'} = '3c2 3c3 2c2 hybrid';
$ref_length[4]=length('CAGAAAGTTGAAGTCGCAAAA');
$R{4}{'ref'}='CAGAAAGTTGAAGTCGCAAAA'; #1c1, 2c1, 3c1

# Region 5
$R{5}{'GGCGTCGTTACCGCC'} = 'ref 1c1 1c2 1c3 2c1 2c2 2c5 2c6 3c1 3c2 3c3 6c2 6c3 7c1 uss'; 
$R{5}{'GGCGTCGTCACCGCC'} = 'v282 1c4 1c5 2c3 2c4 6c1';
$ref_length[5]=length('GGCGTCGTTACCGCC');
$R{5}{'ref'}='GGCGTCGTTACCGCC'; # 1c2, 2c1, 3c1

# Region 6
$R{6}{'GAAATGAAACCAAGCGG'} = '1c1';
$R{6}{'CAAATGAATCCAAGCGG'} = '1c2';
$R{6}{'CAAATGGCTTCAAGCAA'} = '1c5';
$R{6}{'ACAATGCTTTCAAGCGG'} = '1c3 2c2';
$R{6}{'ACAATGAATTCAAGCAA'} = '2c3';
$R{6}{'GAAATGGCTTCAACCGG'} = '1c4 2c4 2c5';
$R{6}{'ACAATGGCTTCAAGCAA'} = '2c6';
$R{6}{'AAAATGCTTTCAAGCGG'} = '3c3 6c3';
$R{6}{'GAAATGAAATCAGACGG'} = '6c2 7c1';
$R{6}{'GAAATGGCTTCAAGCGG'} = 'uss';
$R{6}{'CAAATGGCTTCAACCGG'} = 'ref 2c1 3c1 3c2 6c1'; # 2c1, 3c1, 3c2, 6c1
$R{6}{'CAAATGACTTCAACCGG'} = 'false';
$R{6}{'CAAATGGCTTCAAACCGG'} = 'indel';
$ref_length[6]=length('CAAATGGCTTCAACCGG');
$R{6}{'ref'}='CAAATGGCTTCAACCGG'; # 2c1, 3c1, 3c2, 6c1

# Region 7
$R{7}{'TGAAATCAAAGACAAAAA'} = '1c2';
$R{7}{'AGAAATCAAAGACAAAAA'} = '1c5';
$R{7}{'AGAAATCAAAGACAAAAG'} = '2c3';
$R{7}{'TGAAATCAAAGGCAAAAA'} = '1c3 1c4 2c2 2c4 6c1';
$R{7}{'AGAAATCAAAGGCAAAAA'} = '1c1 3c2 6c2 v322/326';
$R{7}{'AGAAATCCAAGGCAAAAA'} = '2c5 3c1 6c3 vA362G';
$R{7}{'AGAAATCCAAGGCAAAAG'} = '2c6 3c3 7c1';
$R{7}{'AGAAATCAAAGGCAAAAG'} = 'uss';
$R{7}{'AGAAATCCAAGACAAAAA'} = 'ref 2c1'; #2c1 
$R{7}{'AGAAATACCAAGACAAAAA'} = 'indel';
$R{7}{'AGAAATCCAAGCAAAAA'} = 'indel';
$R{7}{'AGAAATCCAAAGACAAAAA'} = 'indel';
$R{7}{'AGAAAATCCAAGACAAAAA'} = 'indel';
$R{7}{'TGAAATCCAAGACAAAAA'} = 'var (A315T)';
$R{7}{'AGAAATCCAAAGGCAAAAGA'} = '2c6 3c3 7c1 hybrid';
$ref_length[7]=length('AGAAATCCAAGACAAAAA');
$R{7}{'ref'}='AGAAATCCAAGACAAAAA'; #2c1 
$R{7}{'GAAATCCAAAGGCAAAAA'}='var';

# Region 8
$R{8}{'CCAAGCGTGAAGACGGTTC'} = '1c1';
$R{8}{'CCAGGCGTGAAAACGGTTC'} = '1c3';
$R{8}{'CCAAGCGTCAAGACGGTC'}  = '1c4';
$R{8}{'CCAAGCGTCAAGCCGGTTC'} = '2c2';
$R{8}{'CCAAGCGTGAAAACGGTTC'} = '1c2 1c5 2c3 2c5 v355/358';
$R{8}{'CCAAGCGTGAAAACGGTC'}  = '2c6';
$R{8}{'GCAGGCGTGAAAACGGTTC'} = '6c2';
$R{8}{'CCAAGCGTGAAGCCGGTTC'} = '3c3 6c3 hybrid 7c1';
$R{8}{'CCAGGCGTGAAGCCGGTTC'} = '7c1';
$R{8}{'CCAGGCGTCAAGACGGTTC'} = '3c2 uss vA350G';
$R{8}{'CCAAGCGTCAAGACGGTTC'} = 'ref 2c1 2c4 3c1 6c1'; # 2c1 2c4 3c1 6c1
$R{8}{'CCAAAGCGTCAAGACGGTTC'} = 'indel';
$ref_length[8]=length('CCAAGCGTCAAGACGGTTC');
$R{8}{'ref'}='CCAAGCGTCAAGACGGTTC'; # 2c1 2c4 3c1 6c1

# Region 9
$R{9}{'AAGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAA'}                   = '1c1';
$R{9}{'ACGCGCAACGACGCCAAAGCCGACGCCAAAGACGACACCGTCACCGCCATCGAA'}                         = '1c2';
$R{9}{'ACGCGCGCCAAAGCCGACGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC'}                   = '1c3';
$R{9}{'AAGCGCGACGCCGGCGCCAAAACCGGCGCCGACGACGTCAAAGCCGACGGCAAAGACACCGACAAAATCAAC'}       = '1c4';
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
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'ref 2c1'; #2c1
$R{9}{'CAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = '2c2';
$R{9}{'AAGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}                         = 'var';
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAA'}                         = 'var';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCGCCAACGCCATCGAC'}                          = 'var';  # added Oct 15, 2013
$ref_length[9]=length('ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC');
$R{9}{'ref'}='ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'; #2c1
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGATGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGTCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGTGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGTGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACTGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAATGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACACCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCAACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ATGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACACCAACAACGCCATCGAC'} = "var";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCATTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAAACAACGCCATCGAC'} = "indel";
$R{9}{'ACGCACACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCAACGACGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAC'} = "3c1 hybrid";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCTGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACTGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'} = "false";
$R{9}{'ACGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATTGAC'} = "false";

$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGTCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGACCGACGCCAACAACGCCATCGAC'}='indel';
$R{9}{'CACGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCGCCAACGCCATCGAA'}='1c1 3c1 hybrid';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGTCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCTATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACATCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCAACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACGACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAAACGCCATCGA'}='indel';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAT'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCTAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGATGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCATCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCAGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGATAACGACGACACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAACGCCATCAAC'}='var';
$R{9}{'CGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCC'}='1c1 hybrid';
$R{9}{'CGCGCACCGGCGACAACGACGATACCGTTGCCGACGCCAACAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAACAATGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAATAACGCCATCGAC'}='mismatch';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTACCGACGCCAACAACGCCATCGAC'}='var';

$R{9}{'AAGCGCGACGCCGGCGCCAAAGCCGACGACGTCAAAGCCGACGCCGCCAACGCCATCGAC'} = '1c1 hybrid'; #see github issue 24
 
# Region 10
# If the regions are shorter than the reference, then they just need a "TA" at the end
$R{10}{'CGATGAATCATCTGCCGTT'}    = '1c3';
$R{10}{'CGATAAATCATCTGCCGTT'}    = '1c4';
$R{10}{'TGATACGTCATCTGCCAAA'}    = '2c2';
$R{10}{'TGATACGTCATCTGCCACCTA'}  = '2c2'; # added Oct 15, 2013
$R{10}{'CGATAAATCATCTGCCACCTA'}  = '2c3 hybrid 7c1'; #also a var
$R{10}{'CGATGAATCATCGTTGCCGG'}   = '1c5 2c5 2c6';
$R{10}{'CGACCCGTTCTCTGCTAGC'}    = '3c3';
$R{10}{'CGATAAACATGATGCCAAA'}    = '2c1 6c1'; #2c1 has "TG" appended
$R{10}{'CGATGAATCATCTGCCGTTTA'}  = '6c2';
$R{10}{'CGATAAATCAACTGCCGTT'}    = '3c2 6c3';
$R{10}{'CGATAAATCAACTGCCAAA'}    = '7c1';
$R{10}{'CGATAAAATCAACTGCCAAA'}   = '7c1';
$R{10}{'CGATGAATCAACTGCCAAA'}    = '7c1'; # added Oct 15, 2013
$R{10}{'CGATGAACCAACTGCCACCTA'}  = 'uss';
$R{10}{'CGATGAATCATCTGCCACCTA'}  = 'ref 1c1 1c2 2c4 3c1'; #1c1 1c2 2c4 3c1
$R{10}{'CGATGAATCATCTGCCAAATA'}  = 'var491';
$R{10}{'CGATGAATCATCTGTCACCTA'}  = 'false';
$R{10}{'CGATGAATCATCTGCCACCCTA'} = 'indel';
$R{10}{'CGATGAAATCATCTGCCACCTA'} = 'indel';
$R{10}{'CGATGAATCATCTGCCCACCTA'} = 'indel';
$R{10}{'CGATAAATCAAACTGCCAAATA'} = '7c1';
$ref_length[10]=length('CGATGAATCATCTGCCACCTA');
$R{10}{'ref'}='CGATGAATCATCTGCCACCTA'; #1c1 1c2 2c4 3c1
$R{10}{'CGATAAATCAACTGCCACCTA'}  = 'var';
