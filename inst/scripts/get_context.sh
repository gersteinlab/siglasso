#!/bin/bash

VCFs="$1"
OUTPUT="$2"
CONTEXT_LEN="$3"
BEDTOOLS="$4"
REF_GENOME="$5"

# check the chromosome format in REF_GENOME
chr_name=`head -1 $REF_GENOME | cut -f1`


while IFS='' read -r line || [[ -n "$line" ]]; do
	VCF=`echo $line | awk -F"[ |\t|,]" '{print $1}'`
	SAMPLE=`echo $line | awk -F"[ |\t|,]" '{if (NF > 1) {print $2} \
			else {gsub(".*\/",""); print $1}}'`
    if [[ $chr_name == *"chr"* ]]
    then
	    grep  -v "#" $VCF | awk -v var=$CONTEXT_LEN \
		    	'{OFS="\t"; print $1, $2-1-var, $2+var}' | \
			    $BEDTOOLS getfasta -fi $REF_GENOME -bed - | \
			    grep -v "^>" | awk '{print toupper($0)}' > $OUTPUT".tmp"
    else
	    grep  -v "#" $VCF | sed "s/^chr//g" | awk -v var=$CONTEXT_LEN \
		    	'{OFS="\t"; print $1, $2-1-var, $2+var}' | \
			    $BEDTOOLS getfasta -fi $REF_GENOME -bed - | \
			    grep -v "^>" | awk '{print toupper($0)}' > $OUTPUT".tmp"
        
    fi 
	grep  -v "#" $VCF | awk -v var=$SAMPLE '{print $5"\t"var}' | \
			paste -d"\t" $OUTPUT".tmp" - >>  $OUTPUT
	rm $OUTPUT".tmp"
done < "$VCFs"
