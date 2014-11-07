#!/bin/bash

# run ./run.sh and make the associated latex file.
# usage: 
# ./run2.sh <input filename> <bool for whether to round or not>

if [ "$1" != "" ]; then

    ./run.sh $1 $2

    echo "\def\filename{$1}" > tex/defrun.tex
    
    cd tex
    echo "running latex..."
    latexmk -pdf transforms &> /dev/null
    echo "   done."
    cd -
    
    outfile=$(echo $1 | sed 's/dinputs//'| sed 's/cinputs//')

    mkdir -p xcorr_pdfs
    cp tex/transforms.pdf xcorr_pdfs/$outfile.pdf


else
    echo "specify the file!"
fi
