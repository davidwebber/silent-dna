#!/bin/bash
today=$(shell date +%b%d)
#./SilentDNA.pl SampleFile.txt     > sample_$today.csv

all:  MID10_FA1090_${today}.csv  MID31_FA1090_${today}.csv    MID46_FA1090_${today}.csv MID1_3_FA1090_${today}.csv   MID33_FA1090_${today}.csv    MID47_FA1090_${today}.csv MID1_4_FA1090_${today}.csv   MID35_FA1090_${today}.csv    MID4_FA1090_${today}.csv MID14_FA1090_${today}.csv    MID36_FA1090_${today}.csv    MID52_FA1090_${today}.csv MID15_FA1090_${today}.csv    MID37_3_FA1090_${today}.csv  MID54_FA1090_${today}.csv MID16_3_FA1090_${today}.csv  MID37_4_FA1090_${today}.csv  MID58_3_FA1090_${today}.csv MID16_4_FA1090_${today}.csv  MID39_FA1090_${today}.csv    MID58_4_FA1090_${today}.csv MID21_FA1090_${today}.csv    MID40_3_FA1090_${today}.csv  MID59_FA1090_${today}.csv MID23_FA1090_${today}.csv    MID40_4_FA1090_${today}.csv  MID62_FA1090_${today}.csv MID24_FA1090_${today}.csv    MID41_FA1090_${today}.csv    MID67_FA1090_${today}.csv MID25_FA1090_${today}.csv    MID42_FA1090_${today}.csv    MID74_FA1090_${today}.csv MID28_FA1090_${today}.csv    MID43_FA1090_${today}.csv MID30_FA1090_${today}.csv    MID45_FA1090_${today}.csv

MID10_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID10_FA1090.fna   > MID10_FA1090_${today}.csv
MID1_3_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID1_3_FA1090.fna  > MID1_3_FA1090_${today}.csv
MID1_4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID1_4_FA1090.fna  > MID1_4_FA1090_${today}.csv
MID14_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID14_FA1090.fna   > MID14_FA1090_${today}.csv
MID15_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID15_FA1090.fna   > MID15_FA1090_${today}.csv
MID16_3_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID16_3_FA1090.fna > MID16_3_FA1090_${today}.csv
MID16_4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID16_4_FA1090.fna > MID16_4_FA1090_${today}.csv
MID21_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID21_FA1090.fna   > MID21_FA1090_${today}.csv
MID23_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID23_FA1090.fna   > MID23_FA1090_${today}.csv
MID24_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID24_FA1090.fna   > MID24_FA1090_${today}.csv
MID25_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID25_FA1090.fna   > MID25_FA1090_${today}.csv
MID28_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID28_FA1090.fna   > MID28_FA1090_${today}.csv
MID30_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID30_FA1090.fna   > MID30_FA1090_${today}.csv
MID31_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID31_FA1090.fna   > MID31_FA1090_${today}.csv
MID33_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID33_FA1090.fna   > MID33_FA1090_${today}.csv
MID35_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID35_FA1090.fna   > MID35_FA1090_${today}.csv
MID36_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID36_FA1090.fna   > MID36_FA1090_${today}.csv
MID37_3_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID37_3_FA1090.fna > MID37_3_FA1090_${today}.csv
MID37_4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID37_4_FA1090.fna > MID37_4_FA1090_${today}.csv
MID39_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID39_FA1090.fna   > MID39_FA1090_${today}.csv
MID40_3_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID40_3_FA1090.fna > MID40_3_FA1090_${today}.csv
MID40_4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID40_4_FA1090.fna > MID40_4_FA1090_${today}.csv
MID41_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID41_FA1090.fna   > MID41_FA1090_${today}.csv
MID42_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID42_FA1090.fna   > MID42_FA1090_${today}.csv
MID43_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID43_FA1090.fna   > MID43_FA1090_${today}.csv
MID45_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID45_FA1090.fna   > MID45_FA1090_${today}.csv
MID46_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID46_FA1090.fna   > MID46_FA1090_${today}.csv
MID47_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID47_FA1090.fna   > MID47_FA1090_${today}.csv
MID4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID4_FA1090.fna    > MID4_FA1090_${today}.csv
MID52_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID52_FA1090.fna   > MID52_FA1090_${today}.csv
MID54_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID54_FA1090.fna   > MID54_FA1090_${today}.csv
MID58_3_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID58_3_FA1090.fna > MID58_3_FA1090_${today}.csv
MID58_4_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID58_4_FA1090.fna > MID58_4_FA1090_${today}.csv
MID59_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID59_FA1090.fna   > MID59_FA1090_${today}.csv
MID62_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID62_FA1090.fna   > MID62_FA1090_${today}.csv
MID67_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID67_FA1090.fna   > MID67_FA1090_${today}.csv
MID74_FA1090_${today}.csv: SilentDNA.pl
	./SilentDNA.pl MID74_FA1090.fna   > MID74_FA1090_${today}.csv
