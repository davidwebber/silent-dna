#!/bin/bash
today=`date +%b%d`
./SilentDNA.pl SampleFile.txt   > sample_$today.csv
./SilentDNA.pl MID4_FA1090.fna  > MID4_FA1090_$today.csv
./SilentDNA.pl MID24_FA1090.fna > MID24_FA1090_$today.csv
