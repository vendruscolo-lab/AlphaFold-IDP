#!/bin/bash
for i in $(echo 5);do
	echo $i > k.dat
	 python pulchra.py 
done
