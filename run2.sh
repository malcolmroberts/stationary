#!/bin/bash

# run ./run.sh and make the associated latex file.
# usage: 
# ./run2.sh <input filename> <bool for whether to round or not>

if [ "$1" != "" ]; then

    ./run.sh $1 $2

    echo "\def\filename{$1}" > tex/defrun.tex
    if [ "$2" != "" ]; then
	echo "\def\aval{$2}" > tex/def_a.tex
    else
	echo "\def\aval{start}" > tex/def_a.tex
    fi
    if [ "$3" != "" ]; then
	echo "\def\bval{$3}" > tex/def_b.tex
    else
	echo "\def\bval{end}" > tex/def_b.tex
    fi
    
    
    cd tex
    latexmk -pdf transforms
    cd -
    
    outfile=$(echo $1 | sed 's/dinputs//'| sed 's/cinputs//')

    echo $outfile
    mkdir -p cinputs/xcorr_pdfs
    cp tex/transforms.pdf cinputs/xcorr_pdfs/$outfile.pdf

    echo $1

else
    echo "specify the file!"
fi
