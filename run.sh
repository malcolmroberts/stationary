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

# original file minus linear term, and add the typical cycle
asy -f pdf plot.asy  -u "filenames=\"data,data.typ\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""
mv plot.pdf data.pdf

# autocorrelation
asy -f pdf plot.asy  -u "filenames=\"data.ac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"linlin\""
mv plot.pdf data_ac.pdf

# FFT of autocorrelation
asy -f pdf plot.asy  -u "filenames=\"data.fac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"loglog\""
mv plot.pdf data_fac.pdf

nperiods=$(cat nperiods)

if [ $nperiods != "1" ]; then
    # find all the typical run files, turn newlines into commas,
    # remove last comma.
    TYPS=$(ls -1 | egrep 'data.ytyp[0-9]' | tr '\n' ','| sed s'/.$//' )
    asy -f pdf plot.asy  -u "filenames=\"${TYPS}\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""
    mv plot.pdf data_typ.pdf
fi


# TODO: compile the pdf

# TODO: put the filename in the PDF

# TODO: add the error term somwehere in the PDF,
