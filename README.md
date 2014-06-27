silent-dna
==========

Which silent copy did this read come from?


##Installation=

###Easy
$ sudo apt-get install bioperl

###Hard
Follow the instructions at 
http://www.bioperl.org/wiki/Installing_BioPerl
under the headings appropriate for your operating system.
(I tried http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix under sections
 - Preparing to Install
 - Installing using CPAN (first do the steps under "Installing in a personal module area")
)

##Running on FA1090
make -j 4 FA1090
./check_completion.sh
./remove_short.sh OR ./remove_short_VD300.sh (edit the shell script to change the day)
./count_sequences.sh OR ./count_sequences_VD300-and-like.sh | tee counted_sequences.txt (edit the shell script to change the day)

optional:
./boolean_sequences.sh
./interpret_boolean.sh

##Getting the tails of sequences
make -f Makefile.tails
cat *VD300-like.tail | sort | uniq -c | sort -n > VD300-like-tails.txt
