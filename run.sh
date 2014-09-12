#!/bin/bash

# A simple bash script to plot 

if [ "$1" == "" ]; then
    echo "Gotta specify a file!"
    exit
fi

echo "Running " $1

./change_format.py $1 cfile

./stationary.py -f cfile

# original file minus linear term, and add the typical cycle
asy -f pdf plot.asy  -u "filenames=\"data,data.typ\"; xlabel=\"time\"; ylabel=\"signal\"; sscale=\"linlin\""
mv plot.pdf data.pdf

# autocorrelation
asy -f pdf plot.asy  -u "filenames=\"data.ac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"linlin\""
mv plot.pdf data_ac.pdf


# FFT of autocorrelation
asy -f pdf plot.asy  -u "filenames=\"data.fac\"; xlabel=\"lag\"; ylabel=\"correlation\"; sscale=\"loglog\""
mv plot.pdf data_fac.pdf


# TODO: compile the pdf

# TODO: put the filename in the PDF

# TODO: add the error term somwehere in the PDF,