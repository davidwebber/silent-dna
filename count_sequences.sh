echo --- Region  1 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[3] }' | sort | uniq -c
echo --- Region  2 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[4] }' | sort | uniq -c
echo --- Region  3 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[5] }' | sort | uniq -c
echo --- Region  4 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[6] }' | sort | uniq -c
echo --- Region  5 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[7] }' | sort | uniq -c
echo --- Region  6 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[8] }' | sort | uniq -c
echo --- Region  7 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[9] }' | sort | uniq -c
echo --- Region  8 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[10] }' | sort | uniq -c
echo --- Region  9 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[11] }' | sort | uniq -c
echo --- Region 10 ---
	cat sample_Nov19.csv | awk '{ split($0,a,","); print a[12] }' | sort | uniq -c
