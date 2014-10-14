#!/bin/bash

# A simple bash script to generate output from stationary.py
# usage:
# ./run.sh <filename> <bool for whether to round or not>

if [ "$1" == "" ]; then
    echo "Gotta specify a file!"
    exit
fi

rstring=""
if [ "$2" != "" ]; then
    rstring="-r "$2
fi

echo "Running " $1

./change_format.py $1 cfile

rm -f data.ytyp*

./stationary.py -f cfile $rstring

startval=$(cat startval)

# original file minus linear term, and add the typical cycle
asy -f pdf plot.asy  -u "filenames=\"data.in\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\" ; start=$startval"
mv plot.pdf data.pdf

if [ "$startval" -ge "0" ]; then

    # autocorrelation
    asy -f pdf plot.asy  -u "filenames=\"data.ac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"linlin\" ; x95=true"
    mv plot.pdf data_ac.pdf

    # FFT of autocorrelation
    asy -f pdf plot.asy  -u "filenames=\"data.fac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"loglog\""
    mv plot.pdf data_fac.pdf

    nperiods=$(cat nperiods)

    if [ $nperiods != "0" ]; then
	# find all the typical run files, turn newlines into commas,
	# remove last comma.
	TYPS=$(ls -1 | egrep 'data.ytyp[0-9]' | tr '\n' ','| sed s'/.$//' )
	asy -f pdf plot.asy  -u "filenames=\"${TYPS}\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""
	mv plot.pdf data_typ.pdf
    fi
fi
