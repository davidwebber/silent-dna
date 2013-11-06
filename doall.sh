#!/bin/bash
today=`date +%b%d`
./SilentDNA.pl MID40_3_FA1090.fna > MID40_3_FA1090_$today.csv
./SilentDNA.pl MID40_4_FA1090.fna > MID40_4_FA1090_$today.csv
