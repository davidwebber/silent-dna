#!/bin/bash
today=`date +%b%d`
./SilentDNA.pl SampleFile.txt     > sample_$today.csv
./SilentDNA.pl MID04_FA1090.fna   > MID04_FA1090_$today.csv
./SilentDNA.pl MID10_FA1090.fna   > MID10_FA1090_$today.csv
./SilentDNA.pl MID24_FA1090.fna   > MID24_FA1090_$today.csv
./SilentDNA.pl MID40_3_FA1090.fna > MID40_3_FA1090_$today.csv
./SilentDNA.pl MID40_4_FA1090.fna > MID40_4_FA1090_$today.csv
