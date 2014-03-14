#for n in `tail *.fna | grep ^\> | awk '{print substr($1,2,30)}' `; do echo $n; grep $n *.csv; done
for n in `tail MID*.fna | grep ^\> | awk '{print substr($1,2,30)}' `; do grep $n *.csv; done | wc
