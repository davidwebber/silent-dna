for n in `ls *.csv`
do
  echo $n
  cat $n | awk '{ split($0,a,","); print a[16] }' | sort | uniq -c | sort -n
  wc $n | awk '{print "   "$1" total"}'
done