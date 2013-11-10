#!/bin/bash
today=`date +%b%d`
./SilentDNA.pl SampleFile.txt   > sample_$today.csv
./SilentDNA.pl MID04_FA1090_shortOnly.fna > MID04_shortOnly_$today.csv
./SilentDNA.pl MID04_FA1090.fna > MID04_FA1090_$today.csv
