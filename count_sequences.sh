day=Mar04
echo --- Region  1 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[3] }'  | sort | uniq -c | sort -n
echo --- Region  2 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[4] }'  | sort | uniq -c | sort -n
echo --- Region  3 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[5] }'  | sort | uniq -c | sort -n
echo --- Region  4 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[6] }'  | sort | uniq -c | sort -n
echo --- Region  5 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[7] }'  | sort | uniq -c | sort -n
echo --- Region  6 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[8] }'  | sort | uniq -c | sort -n
echo --- Region  7 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[9] }'  | sort | uniq -c | sort -n
echo --- Region  8 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[10] }' | sort | uniq -c | sort -n
echo --- Region  9 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[11] }' | sort | uniq -c | sort -n
echo --- Region 10 ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[12] }' | sort | uniq -c | sort -n
echo --- nVar       ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[13] }' | sort | uniq -c | sort -n
echo --- empty      ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[14] }' | sort | uniq -c | sort -n
echo --- verdict    ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[15] }' | sort | uniq -c | sort -n
echo --- altVerdict ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[16] }' | sort | uniq -c | sort -n
echo --- P_minus ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[17] }' | sort | uniq -c | sort -n
echo --- basepairs ---
	cat filtered_*_$day.csv | awk '{ split($0,a,","); print a[18] }' | sort | uniq -c | sort -n
