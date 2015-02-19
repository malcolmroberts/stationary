#!/bin/bash

# Create multiple pdf files based on input filenames.
# usage:
# ./run3.sh <bool for whether to round or not> <stats> <p>

#datadir="dinputs/"
datadir="cinputs/"

rstring=""
if [ "$1" != "" ]; then
    rstring=$1
fi

#INPUTS=$(find $datadir ! -name '*arm*' ! -name '*pdf*' )
# Need two calls to find to get the files in the right order
INPUTS=$(find $datadir -type f -name No[0-9]_* -not -name '*pdf*' | sort  )
INPUTS2=$(find $datadir -type f -name No[0-9][0-9]_* -not -name '*pdf*' | sort  )

INPUTS=${INPUTS}" "
INPUTS=${INPUTS}${INPUTS2}

mkdir -p xcorr_pdfs

rm start_of_stationarity.csv
echo -e "#start_of_stationarity" > start_of_stationarity.csv

rm num_periods.csv
echo -e "#number of periods" > num_periods.csv
rm period_lengths.csv
echo -e "#length of first cycle" > period_lengths.csv
rm period_powers.csv
echo -e "#power of first cycle" > period_powers.csv

for i in $INPUTS:
do
    i=$(echo $i | sed 's/://')
    if [ "$i" != $datadir ]; then
	echo $i
	ionice -c 3 nice -n 19 ./run2.sh $i $rstring $2 $3
	
	startval=$(cat startval)
	echo -e $startval >> start_of_stationarity.csv

	nper=$(cat nperiods)
	echo -e $nper >> num_periods.csv
	rm nperiods

	perlength=$(cat period_length.csv)
	echo -e $perlength >> period_lengths.csv
	rm period_length.csv

	dapowah=$(cat period_power.csv)
	echo -e $dapowah >> period_powers.csv
	rm period_power.csv
    fi
done

asy -f pdf onset.asy
mv onset.pdf xcorr_pdfs/
mv num_periods.csv xcorr_pdfs/
mv start_of_stationarity.csv xcorr_pdfs/
mv period_lengths.csv xcorr_pdfs/
mv period_powers.csv xcorr_pdfs/

