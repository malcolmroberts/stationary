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

rm -f output/start_of_stationarity.csv
touch output/start_of_stationarity.csv
echo -e "#start_of_stationarity" > output/start_of_stationarity.csv

rm -f output/num_periods.csv
touch output/num_periods.csv
echo -e "#number of periods" > output/num_periods.csv

rm -f output/period_lengths.csv
touch output/period_lengths.csv
echo -e "#length of first cycle" > output/period_lengths.csv

rm -f output/period_powers.csv
touch output/period_powers.csv
echo -e "#power of first cycle" > output/period_powers.csv

rm -f output/nonperiod_powers.csv
touch output/nonperiod_powers.csv
echo -e "#power of non-periodic part" > output/nonperiod_powers.csv

for i in $INPUTS:
do
    i=$(echo $i | sed 's/://')
    if [ "$i" != $datadir ]; then
	echo $i
	ionice -c 3 nice -n 19 ./run2.sh $i $rstring $2 $3
	
	startval=$(cat output/startval)
	echo -e $startval >> output/start_of_stationarity.csv
	rm -f output/startval

	nper=$(cat output/nperiods)
	echo -e $nper >> output/num_periods.csv
	rm -f output/nperiods

	perlength=$(cat output/period_length.csv)
	echo -e $perlength >> output/period_lengths.csv
	rm -f output/period_length.csv

	period_power=$(cat output/period_power.csv)
	echo -e $period_power >> output/period_powers.csv
	rm -f output/period_power.csv

	nonperiod_power=$(cat output/nonperiod_power.csv)
	echo -e $nonperiod_power >> output/nonperiod_powers.csv
	rm -f output/nonperiod_power.csv
    fi
done

asy -f pdf onset.asy
mv onset.pdf xcorr_pdfs/

asy -f pdf detected_period_power.asy -u" period_powers_filename = \"output/period_powers.csv\"; nonperiod_powers_filename = \"output/nonperiod_powers.csv\""
mv detected_period_power xcorr_pdfs/

asy -f pdf detected_period_length.asy -u" period_lengths_filename = \"output/period_lengths.csv\""

mv detected_period_length.pdf xcorr_pdfs/

mv num_periods.csv xcorr_pdfs/
mv output/start_of_stationarity.csv xcorr_pdfs/
mv output/period_lengths.csv xcorr_pdfs/
mv output/period_powers.csv xcorr_pdfs/
mv output/nonperiod_powers.csv xcorr_pdfs/

