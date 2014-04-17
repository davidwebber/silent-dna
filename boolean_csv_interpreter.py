#!/usr/bin/env python

import sys
import csv

region_length=[
    181-159+1, #make the computer do the subtraction
    210-187+1,
    239-223+4,
    273-253+1,    
    288-274+1,
    305-289+1,
    332-315+1,
    365-347+1,
    447-394+25,
    494-474+1];

#print region_length;

if __name__ == "__main__":
	#print "hello world"
	if len(sys.argv)<2:
		print "Usage:", sys.argv[0], "<input filename>"
		sys.exit(0)

	with open(sys.argv[1]) as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			tot=0
			for i in range(0,len(row)):
				if (row[i]=='1'):
					print i+1,
					tot += region_length[i]
			print ","+str(tot)

			#for (i in range (0,len(row))):
			#	print i, "mississippi"