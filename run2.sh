#!/bin/bash

# run ./run.sh and make the associated latex file.
# usage: 
# ./run2.sh <input filename> <bool for whether to round or not>

if [ "$1" != "" ]; then

    ./run.sh $1 $2

    echo "\def\filename{$1}" > tex/defrun.tex
    
    cd tex
    latexmk -pdf transforms &> /dev/null
    cd -
    
    outfile=$(echo $1 | sed 's/dinputs//'| sed 's/cinputs//')

    echo $outfile
    mkdir -p cinputs/xcorr_pdfs
    cp tex/transforms.pdf cinputs/xcorr_pdfs/$outfile.pdf

    echo $1

else
    echo "specify the file!"
fi
