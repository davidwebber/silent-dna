#!/bin/bash

day=May27

cat MID11_VD300_$day.csv        | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID11_VD300_$day.csv  
cat MID13_VD300-like_$day.csv   | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID13_VD300-like_$day.csv 
cat MID17_VD300-like_$day.csv   | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID17_VD300-like_$day.csv 
cat MID19_VD300-like_$day.csv   | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID19_VD300-like_$day.csv  
cat MID20_3_VD300-like_$day.csv | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID20_3_VD300-like_$day.csv
cat MID20_4_VD300-like_$day.csv | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID20_4_VD300-like_$day.csv
cat MID32_VD300-like_$day.csv   | awk '{split($0,a,","); if (a[18]>=480){print $0}}' > filtered_MID32_VD300-like_$day.csv  
