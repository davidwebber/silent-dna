day=Mar04
echo --- Region  1 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[3] }'  | sort | uniq -c | sort -n
echo --- Region  2 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[4] }'  | sort | uniq -c | sort -n
echo --- Region  3 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[5] }'  | sort | uniq -c | sort -n
echo --- Region  4 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[6] }'  | sort | uniq -c | sort -n
echo --- Region  5 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[7] }'  | sort | uniq -c | sort -n
echo --- Region  6 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[8] }'  | sort | uniq -c | sort -n
echo --- Region  7 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[9] }'  | sort | uniq -c | sort -n
echo --- Region  8 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[10] }' | sort | uniq -c | sort -n
echo --- Region  9 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[11] }' | sort | uniq -c | sort -n
echo --- Region 10 ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[12] }' | sort | uniq -c | sort -n
echo --- nVar       ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[13] }' | sort | uniq -c | sort -n
echo --- empty      ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[14] }' | sort | uniq -c | sort -n
echo --- verdict    ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[15] }' | sort | uniq -c | sort -n
echo --- altVerdict ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[16] }' | sort | uniq -c | sort -n
echo --- P_minus ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[17] }' | sort | uniq -c | sort -n
echo --- basepairs ---
	cat short_*_$day.csv | awk '{ split($0,a,","); print a[18] }' | sort | uniq -c | sort -n
