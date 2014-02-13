for n in `ls *.csv | colrm 20`; do cp -v ../$n*.csv .; done
