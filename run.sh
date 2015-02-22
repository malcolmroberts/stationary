#!/bin/bash

# A simple bash script to generate output from stationary.py
# usage:
# ./run.sh <filename> <bool for whether to round or not> <stest> <p-value>

if [ "$1" == "" ]; then
    echo "Gotta specify a file!"
    exit
fi

rstring=""
if [ "$2" != "" ]; then
    rstring="-r "$2
fi

sstring=""
if [ "$3" != "" ]; then
    sstring="-s "$3
fi

pstring=""
if [ "$4" != "" ]; then
    pstring="-p "$4
fi

echo "Running " $1

#./change_format.py $1 cfile

rm -f output/*

echo ./stationary.py -f $1 $rstring $sstring $pstring -t 1
./stationary.py -f $1 $rstring $sstring $pstring -t 1

set -x
startval=$(cat output/startval)

# original file minus linear term, and add the typical cycle
asy -f pdf plot.asy  -u "filenames=\"output/data.in\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\" ; start=$startval"
mv plot.pdf data.pdf


if [ "$startval" -ge "0" ]; then

    # autocorrelation
    asy -f pdf plot.asy  -u "filenames=\"output/data.ac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"linlin\" ; x95=true"
    mv plot.pdf data_ac.pdf

    # FFT of autocorrelation
    asy -f pdf plot.asy  -u "filenames=\"output/data.fac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"loglog\""
    mv plot.pdf data_fac.pdf

    nperiods=$(cat output/nperiods)

    if [ $nperiods != "0" ]; then
	# find all the typical run files, turn newlines into commas,
	# remove last comma.
	TYPS=$(ls -1 output/data.ytyp* | egrep 'data.ytyp' | tr '\n' ','| sed s'/.$//' )
	asy -f pdf plot.asy  -u "filenames=\"${TYPS}\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""  &> /dev/null
	mv plot.pdf data_typ.pdf


    fi

    asy -f pdf plot.asy  -u "filenames=\"output/data.np\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""
    mv plot.pdf data_np.pdf

fi
