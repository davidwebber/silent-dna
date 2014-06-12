#This sequence data for VD300-like = ref
#MID13_VD300-like.fna
#MID17_VD300-like.fna
#MID19_VD300-like.fna
#MID20_3_VD300-like.fna
#MID20_4_VD300-like.fna
#MID32_VD300-like.fna

#VD300 keeps its tag.  VD300-like gets tagged "ref"

# offset between reference sequence and reads

$offset=91; 

$nRegions=11;

@region_min=(0, 
    162,
    173,
    195,
    223,
    256,
    285,
    318,
    350,
    398,
    483,
    529);

@region_max=(0,
    168,
    188,
    210,
    242,
    276,
    308,
    335,
    362,
    450,
    507,
    541);

@region_length=(0, 
    $region_max[1]-$region_min[1]+1, #make the computer do the subtraction
    $region_max[2]-$region_min[2]+1,
    $region_max[3]-$region_min[3]+4,
    $region_max[4]-$region_min[4]+1,    
    $region_max[5]-$region_min[5]+1,
    $region_max[6]-$region_min[6]+1,
    $region_max[7]-$region_min[7]+1,
    $region_max[8]-$region_min[8]+1,
    $region_max[9]-$region_min[9]+25,
    $region_max[10]-$region_min[10]+1,
    $region_max[11]-$region_min[11]+1);

$bound[1]  = "CAGCCGT";          #boundary between retions 0 and 1
$bound[2]  = "TATT";             #boundary between regions 1 and 2, original
$bound[3]  = "ATGGCC";           #boundary between regions 2 and 3
$bound[4]  = "GCCGGCGGGCA";      #boundary between regions 3 and 4
$bound[5]  = "AGGCAAATATGTT";    #boundary between regions 3 and 4
$bound[6]  = "GGCGTCGT";         #boundary between regions 3 and 4
$bound[7]  = "CGTAAACAA";        #boundary between regions 6 and 7
$bound[8]  = "ACTCTCCCTGTGGG";   #boundary between regions 7 and 8
#R8- CGGTTCGGTA ... CAGCCGGTTA-R9
$bound[9]  = "CAGCCGGTTA";                 #boundary between regions 8 and 9
#AAACGGTTCG at the end of R8 and CAGCCGGTTA right before R9
$bound[10] = "ACCAAGCACCTGCCGTCAACCTGCCG";    #boundary between regions 9 and 10 for either 1c5 or ref
$bound[11] = "TTAAATTT";

$bitmap{'1c1'}     = 0b0000000000000000000001;
$bitmap{'1c2'}     = 0b0000000000000000000010;
$bitmap{'1c3'}     = 0b0000000000000000000100;
$bitmap{'1c4'}     = 0b0000000000000000001000;
$bitmap{'1c5'}     = 0b0000000000000000010000;
$bitmap{'2c1'}     = 0b0000000000000000100000;
$bitmap{'2c2'}     = 0b0000000000000001000000;
$bitmap{'2c3'}     = 0b0000000000000010000000;
$bitmap{'2c4'}     = 0b0000000000000100000000;
$bitmap{'2c5'}     = 0b0000000000001000000000;
$bitmap{'2c6'}     = 0b0000000000010000000000;
$bitmap{'3c1'}     = 0b0000000000100000000000;
$bitmap{'3c2'}     = 0b0000000001000000000000;
$bitmap{'3c3'}     = 0b0000000010000000000000;
$bitmap{'6c1'}     = 0b0000000100000000000000;
$bitmap{'6c2'}     = 0b0000001000000000000000;
$bitmap{'6c3'}     = 0b0000010000000000000000;
$bitmap{'7c1'}     = 0b0000100000000000000000;
$bitmap{'uss'}     = 0b0001000000000000000000;
$bitmap{'uss2'}    = 0b0010000000000000000000;
$bitmap{'uss3'}    = 0b0100000000000000000000;
$bitmap{'VD300'}   = 0b1000000000000000000000;
$bitmap{'ref'}     = 0b1111111111111111111111;


while (($key, $value) = each %bitmap){
    $reversebitmap{$value} = $key;
}
$reversebitmap{0b0000000000000000000} = 'double_cross';
$reversebitmap{0b0000100000000100000} = '2c1_6c1';

# Region 1 
$R{1}{'CACCGAG'}='ref 2c1 5c1 7c1';
$R{1}{'TGCCGGG'}='1c1 1c3 1c5 uss3';
$R{1}{'TACCGAG'}='1c2 2c2 uss2';
#$R{1}{'TGCCGGG'}='1c3';
$R{1}{'TACCGGG'}='1c4';
#$R{1}{'TGCCGGG'}='1c5';
#1c6
#$R{1}{'CACCGAG'}='2c1';
#$R{1}{'TACCGAG'}='2c2';
#$R{1}{'CACCGAG'}='5c1';
$R{1}{'CACCGAA'}='6c1';
#6c2
#6c3
#$R{1}{'CACCGAG'}='7c1';
#$R{1}{'TACCGAG'}='uss2';
#$R{1}{'TGCCGGG'}='uss3';
$R{1}{'ref'}='CACCGAG'; 

# Region 2
$R{2}{'ACCTGAATCACGGCAA'}='VD300 6c3';
$R{2}{'ACCTGAATCACGGCGA'}='1c1 1c2 2c1 5c1 ref';
#$R{2}{'ACCTGAATCACGGCGA'}='1c2';
$R{2}{'GCCCGAATCACGGCAC'}='1c3 1c5';
$R{2}{'GCCCGAATCACGGCAA'}='1c4';
#$R{2}{'GCCCGAATCACGGCAC'}='1c5';
#1c6
#$R{2}{'ACCTGAATCACGGCGA'}='2c1';
$R{2}{'ACCTGAATCACGGCAT'}='2c2 7c1';
#$R{2}{'ACCTGAATCACGGCGA'}='5c1';
$R{2}{'ACCCGAATAACGGCAA'}='6c1';
#6c2
#$R{2}{'ACCTGAATCACGGCAA'}='6c3';
#$R{2}{'ACCTGAATCACGGCAT'}='7c1';
$R{2}{'ACCTGAATAACGGCGA'}='uss2 uss3';
#$R{2}{'ACCTGAATAACGGCGA'}='uss3';
$R{2}{'ref'}='ACCTGAATCACGGCGA'; 

# Region 3
$R{3}{'GGAAAACAACACTTCT'}='VD300';
$R{3}{'CAAAGACAACGACTCT'}='1c1 5c1';
$R{3}{'GAAAGACAACACTTCT'}='1c2 1c4 2c2 7c1';
$R{3}{'GAAAGACAACGGTGA'}='1c3';
#$R{3}{'GAAAGACAACACTTCT'}='1c4';
$R{3}{'GGAAAACAACGCTTCT'}='1c5 ref';
#1c6
$R{3}{'GGAAGACAACACTTCT'}='2c1';
#$R{3}{'GAAAGACAACACTTCT'}='2c2';
#$R{3}{'CAAAGACAACGACTCT'}='5c1';
$R{3}{'CGCCGACAACGGCGCT'}='6c1';
#6c2
$R{3}{'GGAAAACAAGCCT'}='6c3';
#$R{3}{'GAAAGACAACACTTCT'}='7c1';
$R{3}{'GGAAGACAACGACAAG'}='uss2';
$R{3}{'GGAAGACAACGGCGCT'}='uss3';
$R{3}{'ref'}='GGAAAACAACGCTTCT'; 

# Region 4
$R{4}{'TCCCCCCCCACCGACATCAA'}='VD300';
$R{4}{'TCCCCCCCCACGACATCAAA'}='indel';
$R{4}{'TCCGCTTCAAAAATCAT'}='1c1 5c1 uss3';
$R{4}{'TCCTCCGACAAAATCAA'}='1c2';
$R{4}{'TCCCCCGCCGACAAAATCAA'}='1c3';
$R{4}{'TCCCCCGCCGAAATCAA'}='1c4';
$R{4}{'TCCTCCCCCACCGACATCAA'}='1c5 2c1 ref';
#1c6
#$R{4}{'TCCTCCCCCACCGACATCAA'}='2c1';
$R{4}{'TCTTCTTCATCAATCAA'}='2c2 6c1';
#$R{4}{'TCCGCTTCAAAAATCAT'}='5c1';
#$R{4}{'TCTTCTTCATCAATCAA'}='6c1';
#6c2
$R{4}{'TCCCCCGCCTCCGACATCAA'}='6c3';
$R{4}{'TCCCCCCCCTCCGACATCAA'}='7c1';
$R{4}{'GCCGCCTCCGACATCAA'}='uss2';
#$R{4}{'TCCGCTTCAAAAATCAT'}='uss3';
$R{4}{'TCCGCTTCAACAATCAAAGGCAA'}='uss3 hybrid';
$R{4}{'ref'}='TCCTCCCCCACCGACATCAA'; #1c1, 2c1, 3c1

# Region 5
$R{5}{'AAAGAGGTTGAAGTTAAAAAC'}='VD300';
$R{5}{'AAAGAGGTTGAAGTTAAAAC'}='indel';
$R{5}{'AAAAGCGTTACGGTCGCAAAA'}='ref 1c4 1c5';
$R{5}{'AAGCAAGTTGAAGTCAAAAAC'}='1c1 5c1';
$R{5}{'CAGAAAGTTGAAGTCGCAAAA'}='1c2 1c3';
#$R{5}{'CAGAAAGTTGAAGTCGCAAAA'}='1c3';
#$R{5}{'AAAAGCGTTACGGTCGCAAAA'}='1c4';
#$R{5}{'AAAAGCGTTACGGTCGCAAAA'}='1c5';
#1c6
$R{5}{'CAAAGCGTTACGGTCGCAAAC'}='2c1 6c3 7c1';
$R{5}{'AAGGAAGTTAAAGTCGAAAAC'}='2c2 6c1 uss3';
#$R{5}{'AAGCAAGTTGAAGTCAAAAAC'}='5c1';
#$R{5}{'AAGGAAGTTAAAGTCGAAAAC'}='6c1';
#6c2
#$R{5}{'CAAAGCGTTACGGTCGCAAAC'}='6c3';
#$R{5}{'CAAAGCGTTACGGTCGCAAAC'}='7c1';
$R{5}{'AAAAGCGTTACGGTCGCAAAC'}='uss2';
#$R{5}{'AAGGAAGTTAAAGTCGAAAAC'}='uss3';
$R{5}{'ref'}='AAAAGCGTTACGGTCGCAAAA'; # 1c2, 2c1, 3c1

# Region 6
$R{6}{'TACCGCCACAATGCTTTCAAGCGG'}='VD300 1c6 uss2';
$R{6}{'TACCGCCGAAATGGCTTCAACCGG'}='ref 1c1 1c3 1c5 2c1';
#$R{6}{'TACCGCCGAAATGGCTTCAACCGG'}='1c1';
$R{6}{'TACCGCCACAATGGCTTCAAGCAA'}='1c2';
#$R{6}{'TACCGCCGAAATGGCTTCAACCGG'}='1c3';
$R{6}{'TACCGCCCAAATGAATCCAAGCGG'}='1c4';
#$R{6}{'TACCGCCGAAATGGCTTCAACCGG'}='1c5';
#$R{6}{'TACCGCCACAATGCTTTCAAGCGG'}='1c6';
#$R{6}{'TACCGCCGAAATGGCTTCAACCGG'}='2c1';
$R{6}{'CACCGCCACAATGAATTCAAGCGG'}='2c2';
$R{6}{'CACCGCCACAATGAATTCAAGCAA'}='6c1';
#$R{6}{'TACCGCCCAAATGAAATCAGACGG'}='5c1';
#$R{6}{'CACCGCCACAATGAATTCAAGCGG'}='6c1';
#6c2
$R{6}{'TACCGCCCAAATGAAATCAGACGG'}='6c3 7c1 5c1';
#$R{6}{'TACCGCCCAAATGAAATCAGACGG'}='7c1';
#$R{6}{'TACCGCCACAATGCTTTCAAGCGG'}='uss2';
$R{6}{'CACCGCCCAAATGGCTTCAAGCAA'}='uss3';
$R{6}{'ref'}='TACCGCCGAAATGGCTTCAACCGG'; # 2c1, 3c1, 3c2, 6c1

# Region 7
$R{7}{'TGAAATCAAAGGCAAAAA'}='VD300 uss2';
$R{7}{'AGAAATCCAAGGCAAAAG'}='ref 2c1 2c2';
$R{7}{'AGAAATCCAAGACAAAAACTC'}='var';
$R{7}{'AGAAATCAAAGACAAAAA'}='1c1 1c2 1c3';
#$R{7}{'AGAAATCAAAGACAAAAA'}='1c2';
#$R{7}{'AGAAATCAAAGACAAAAA'}='1c3';
$R{7}{'TGAAATCAAAGACAAAAA'}='1c4';
$R{7}{'AGAAATCAAAGGCAAAAA'}='1c5 1c6 uss3';
#$R{7}{'AGAAATCAAAGGCAAAAA'}='1c6';
#$R{7}{'AGAAATCCAAGGCAAAAG'}='2c1';
#$R{7}{'AGAAATCCAAGGCAAAAG'}='2c2';
$R{7}{'AGAAATCAAAAACAAAAA'}='5c1 6c3 7c1';
$R{7}{'AGAAATCAAAGACAAAAG'}='6c1';
#6c2
#$R{7}{'AGAAATCAAAAACAAAAA'}='6c3';
#$R{7}{'AGAAATCAAAAACAAAAA'}='7c1';
#$R{7}{'TGAAATCAAAGGCAAAAA'}='uss2';
#$R{7}{'AGAAATCAAAGGCAAAAA'}='uss3';
$R{7}{'ref'}='AGAAATCCAAGGCAAAAG'; #2c1

# Region 8
$R{8}{'GCAGGCGTGAAAA'}='VD300 2c2 6c1';
$R{8}{'CCAAGCGTGAAAA'}='ref 1c4 1c5 2c1 5c1 7c1';
$R{8}{'CCAAGCGTCAAGA'}='1c1 1c3 1c6 uss2';
$R{8}{'CCAGGCGTCAAGA'}='1c2';
#$R{8}{'CCAAGCGTCAAGA'}='1c3';
#$R{8}{'CCAAGCGTGAAAA'}='1c4';
#$R{8}{'CCAAGCGTGAAAA'}='1c5';
#$R{8}{'CCAAGCGTCAAGA'}='1c6';
#$R{8}{'CCAAGCGTGAAAA'}='2c1';
#$R{8}{'GCAGGCGTGAAAA'}='2c2';
#$R{8}{'CCAAGCGTGAAAA'}='5c1';
#$R{8}{'GCAGGCGTGAAAA'}='6c1';
#6c2
$R{8}{'CCAGGCGTGAAGC'}='6c3';
#$R{8}{'CCAAGCGTGAAAA'}='7c1';
#$R{8}{'CCAAGCGTCAAGA'}='uss2';
$R{8}{'GCAGGCGTGAAGA'}='uss3';
$R{8}{'ref'}='CCAAGCGTGAAAA'; 

# Region 9
$R{9}{'CGCGCGCCGACGACGACACCGTTGCCGACGCCAAAGACGGCAAAGAAATCGAC'}='VD300';
$R{9}{'CGCGCGGCGCCGGCAACGCCGGCAAAGCCGACGACGTCACCAAAGCCGGCAACGACAACGAAAAAATCAAC'}='1c1 1c3';
$R{9}{'CGCGCGCCAAAGCCAAAGACGCCGACGACGTTACCGACGACGCCGGCACCGACAACGGCGGCAAAGGCAAAATCGAC'}='1c2';
#$R{9}{'CGCGCGGCGCCGGCAACGCCGGCAAAGCCGACGACGTCACCAAAGCCGGCAACGACAACGAAAAAATCAAC'}='1c3';
$R{9}{'CGCGCACCGACGACGCCGCCAAAGACGCCGTTACCGCCGACGCCAAAGACGCCATCGAA'}='1c4';
$R{9}{'CGCGCGCCAAAGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC'}='1c5';
#$R{9}{'CGCGCGCCGCCAAAGACGACGACACCGTTGCCGACGCCAAAGACGGCAAAGAATCGACACCAAGCACCTGCCGTCAACCT'}='5c1 hybrid';
$R{9}{'CGCGCGCCGCCAAAGACGACGACACCGTTGCCGACGCCAAAGACGGCAAAGAATCGAC'}='5c1 hybrid';
$R{9}{'AGCGCAACGACAACGCCGACAACGACGACGTTACCCGCGACGGCACCGACGGCAAAGACAAAATCGAA'}='1c6';
$R{9}{'CGCGCGCCGCCAAAGACGACGACGACGCCGTCACCGCCCGACGGCAACAACAAAATCGAC'}='2c1';
$R{9}{'CGCGCGCCAAAGCCGACGCCGACGCCGACGCCGCCGGCAAAGACACCACCAACATCGAC'}='2c2';
$R{9}{'CGCGCAACGCCAAAGCCAACGACACCGTTGCCGCCGACGGCACCGGCAACGACAAAATCGAA'}='5c1';
$R{9}{'CGCGCAACGCCAACGACGACACCGTCACCGCCGACGGCACCGGCAACGACGGCAAAATCGAC'}='6c1';
$R{9}{'CGCGCGACAACGCCGGCACCGACGCCGTCACCGCCGACACCACCGGCAAAGACAAAGAAATCGAC'}='6c2';
$R{9}{'CGCGCGACAAAGCCGTCACCGACGACGCCGTCAAAGACGTCACCGGCAACGACAAAATCGAA'}='6c3';
$R{9}{'CGCGCGCCGCCAAAGACGACGACGCCGTCACCGCCGACGGCAACAACAAAATCGAC'}='7c1';
$R{9}{'CGCGCACCGGCGACAACGACGACACCGTTGCCGACGCCAAAGACGGCAAAGAAATCGAC'}='uss2';
$R{9}{'AGCGCGACGACACCGCCGCCAAAGCCGGCACCGACGACGTCGCCAAAGACGACGCCGCCGGCAAAAAAATCGAC'}='uss3';
#$R{9}{'CGCGCGCCGCCAAAGACGACGACACGGTTGCCGACGCCAAAGACGGCAAAGAAATCGAC'}='ref';
$R{9}{'CGCGCGCCGCCAAAGACGACGACACCGTTGCCGACGCCAAAGACGGCAAAGAAATCGAC'}='ref';
$R{9}{'ref'}='CGCGCGCCGCCAAAGACGACGACACGGTTGCCGACGCCAAAGACGGCAAAGAAATCGAC'; #2c1

# Region 10
# If the regions are shorter than the reference, then they just need a "TA" at the end
#$R{10}{'CGATAAGGCATCTGATGCCAAATGA'}='VD300';
$R{10}{'CGATAAGGCATCTGATGCCAAATGA'}='ref VD300';
$R{10}{'CGATAACTTTGATGCCAGCTGA'}='1c1';
$R{10}{'CGATAAATCAACTGCCGTT'}='1c2';
$R{10}{'CGATAAATCAACTGCCATT'}='1c3';
$R{10}{'CGATGAATCATCTGCCGGT'}='1c4';
$R{10}{'CGATGAATCATCTGCCGTT'}='1c5';
$R{10}{'CGATGAATCAACTGCCGTTTGCACG'}='1c6';
$R{10}{'CGACACTTCATCTGCCGGTAAGTGA'}='2c1 7c1 uss2';
$R{10}{'CGACGCAGCATCTGCCGTTTGCATA'}='2c2';
$R{10}{'CGATAACTTTGATGCCAGCTG'}='5c1';
$R{10}{'CGACACTTCATCAGCCGGTAAG'}='6c1';
$R{10}{'CGATAAATCATCTGCCGAA'}='6c2';
$R{10}{'CGATGAATCATCTGCCACCTAAGGCAAAT'}='1c4 1c5 uss3 hybrid';
#$R{10}{'CGACACTTCATCTGCCGGTAAGTGA'}='7c1';
#$R{10}{'CGACACTTCATCTGCCGGTAAGTGA'}='uss2';
#$R{10}{'CGAT'}='uss3';
$R{10}{'CGATGAACCAACTGCCACCTAAGGCAAAT'}='uss3';
#$R{10}{'T'}='6c3';
$R{10}{'ref'}='CGATAAGGCATCTGATGCCAAATGA'; #1c11c2 2c4 3c1

$R{11}{'TAAATAAATCAAG'}='ref';  #(both VD300 and VD300like)
$R{11}{'CAAATAAATCAAG'}='5c1 uss3';
$R{11}{'CAAATAAATCAAA'}='1c1';

@ref_length=(0,0,0,0,0,0,0,0,0,0,0);
for ($i=1;$i<=10;$i++){
  $ref_length[$i]=length($R{$i}{'ref'});
}


1; #required for module to load
