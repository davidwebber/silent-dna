=silent-dna=
==========

Which silent copy did this read come from?

==========

=Installation=

==Easy:==
$ sudo apt-get install bioperl

==Hard:==
Follow the instructions at 
http://www.bioperl.org/wiki/Installing_BioPerl
under the headings appropriate for your operating system.
(I tried http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix under sections
 - Preparing to Install
 - Installing using CPAN (first do the steps under "Installing in a personal module area")
)

=Running=
make -j 4
./check_completion.sh
./remove_short.sh OR ./remove_short_VD300.sh
./count_sequences.sh | tee counted_sequences.txt
./boolean_sequences.sh
./interpret_boolean.sh
