#!/bin/bash
cd Dataset/ncbi-genomes-2019-02-19
for f in *.fna.gz; do
	#echo $f;
	name=$(echo "$f" | cut -f 1 -d '.')
	echo $name;
    # do some stuff here with "$f"
    # remember to quote it or spaces may misbehave
done
