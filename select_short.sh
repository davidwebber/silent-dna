#!/bin/bash

day=Mar04

cat MID10_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID10_FA1090_$day.csv  
cat MID1_3_FA1090_$day.csv  | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID1_3_FA1090_$day.csv 
cat MID1_4_FA1090_$day.csv  | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID1_4_FA1090_$day.csv 
cat MID14_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID14_FA1090_$day.csv  
cat MID15_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID15_FA1090_$day.csv  
cat MID16_3_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID16_3_FA1090_$day.csv
cat MID16_4_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID16_4_FA1090_$day.csv
cat MID21_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID21_FA1090_$day.csv  
cat MID23_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID23_FA1090_$day.csv  
cat MID24_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID24_FA1090_$day.csv  
cat MID25_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID25_FA1090_$day.csv  
cat MID28_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID28_FA1090_$day.csv  
cat MID30_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID30_FA1090_$day.csv  
cat MID31_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID31_FA1090_$day.csv  
cat MID33_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID33_FA1090_$day.csv  
cat MID35_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID35_FA1090_$day.csv  
cat MID36_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID36_FA1090_$day.csv  
cat MID37_3_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID37_3_FA1090_$day.csv
cat MID37_4_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID37_4_FA1090_$day.csv
cat MID39_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID39_FA1090_$day.csv  
cat MID40_3_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID40_3_FA1090_$day.csv
cat MID40_4_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID40_4_FA1090_$day.csv
cat MID41_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID41_FA1090_$day.csv  
cat MID42_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID42_FA1090_$day.csv  
cat MID43_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID43_FA1090_$day.csv  
cat MID45_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID45_FA1090_$day.csv  
cat MID46_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID46_FA1090_$day.csv  
cat MID47_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID47_FA1090_$day.csv  
cat MID4_FA1090_$day.csv    | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID4_FA1090_$day.csv   
cat MID52_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID52_FA1090_$day.csv  
cat MID54_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID54_FA1090_$day.csv  
cat MID58_3_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID58_3_FA1090_$day.csv
cat MID58_4_FA1090_$day.csv | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID58_4_FA1090_$day.csv
cat MID59_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID59_FA1090_$day.csv  
cat MID62_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID62_FA1090_$day.csv  
cat MID67_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID67_FA1090_$day.csv  
cat MID74_FA1090_$day.csv   | awk '{split($0,a,","); if (a[18]<480){print $0}}' > filtered_MID74_FA1090_$day.csv  
