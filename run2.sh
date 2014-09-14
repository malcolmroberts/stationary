#!/bin/bash

a_arg=""
if [ "$2" != "" ]; then
    a_arg="$2"
fi
b_arg=""
if [ "$3" != "" ]; then
    b_arg="$3"
fi

if [ "$1" != "" ]; then

    ./run.sh $1 $a_arg $b_arg

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
    
    nperiods=$(cat nperiods)
    
    cd tex
    if [ $nperiods != "1" ]; then
	latexmk -pdf transforms
    else
	latexmk -pdf transforms_np
	mv transforms_np.pdf transforms.pdf
    fi
    cd -
    
    outfile=$(echo $1 | sed 's/dinputs//'| sed 's/cinputs//')

    echo $outfile
    mkdir -p cinputs/xcorr_pdfs
    cp tex/transforms.pdf cinputs/xcorr_pdfs/$outfile.pdf

    echo $1

else
    echo "specify the file!"
fi
