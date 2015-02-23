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

mkdir -p output

rm -f xcorr_pdfs/start_of_stationarity.csv
#touch xcorr_pdfs/start_of_stationarity.csv
echo -e "#start_of_stationarity" > xcorr_pdfs/start_of_stationarity.csv

rm -f xcorr_pdfs/num_periods.csv
touch xcorr_pdfs/num_periods.csv
echo -e "#number of periods" > xcorr_pdfs/num_periods.csv

rm -f xcorr_pdfs/period_lengths.csv
touch xcorr_pdfs/period_lengths.csv
echo -e "#length of first cycle" > xcorr_pdfs/period_lengths.csv

rm -f xcorr_pdfs/period_powers.csv
touch xcorr_pdfs/period_powers.csv
echo -e "#power of first cycle" > xcorr_pdfs/period_powers.csv

rm -f xcorr_pdfs/nonperiod_powers.csv
touch xcorr_pdfs/nonperiod_powers.csv
echo -e "#power of non-periodic part" > xcorr_pdfs/nonperiod_powers.csv

for i in $INPUTS:
do
    i=$(echo $i | sed 's/://')
    if [ "$i" != $datadir ]; then
	echo $i
	ionice -c 3 nice -n 19 ./run2.sh $i $rstring $2 $3
	
	startval=$(cat output/startval)
	echo -e $startval >> xcorr_pdfs/start_of_stationarity.csv
	rm -f output/startval

	nper=$(cat output/nperiods)
	echo -e $nper >> xcorr_pdfs/num_periods.csv
	rm -f output/nperiods

	perlength=$(cat output/period_length.csv)
	echo -e $perlength >> xcorr_pdfs/period_lengths.csv
	rm -f output/period_length.csv

	period_power=$(cat output/period_power.csv)
	echo -e $period_power >> xcorr_pdfs/period_powers.csv
	rm -f output/period_power.csv

	nonperiod_power=$(cat output/nonperiod_power.csv)
	echo -e $nonperiod_power >> xcorr_pdfs/nonperiod_powers.csv
	rm -f output/nonperiod_power.csv
    fi
done


asy -f pdf onset.asy  -u" sobs_filename = \"xcorr_pdfs/start_of_stationarity.csv\"; pobs_filename = \"xcorr_pdfs/num_periods.csv\""
mv onset.pdf xcorr_pdfs/

asy -f pdf detected_period_power.asy -u" period_powers_filename = \"xcorr_pdfs/period_powers.csv\"; nonperiod_powers_filename = \"xcorr_pdfs/nonperiod_powers.csv\""
mv detected_period_power.pdf xcorr_pdfs/

asy -f pdf detected_period_length.asy -u" period_lengths_filename = \"xcorr_pdfs/period_lengths.csv\""

mv detected_period_length.pdf xcorr_pdfs/



