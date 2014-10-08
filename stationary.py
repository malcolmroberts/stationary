#!/usr/bin/python -u

import csv # to load the input file.
import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

from cycle import *
from utils import *
from stats import *


# Main program
def main(argv):
    usage="./stationary -f <filename>"

    # Filename for the input data (consisting of tab-separated
    # (time,value) pairs.
    filename=""

    round=False

    # Load the command-line arguments
    try:
        opts, args = getopt.getopt(argv,"f:r:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filename=arg
        if opt in ("-r"):
            round=(arg == "True" or arg == "true" or arg == "1")

    # Check that the file exists (and is not a directory):
    if not (os.path.isfile(filename)):
        print "Error: the specified file \""+filename+ " \" does not exist"
        print usage
        sys.exit(2)
    
    # Read the data
    data = []
    csvReader = csv.reader(open(filename, 'rb'), delimiter='\t')
    for row in csvReader:
        data.append(row)
    print str(len(data))+" data points found"

    # Consider only the stationary part of the data:
    data=stationary_part(data)

    # Put the y-values from the input data into y:
    y=y_part(data)

    # Remove the linear fit:
    ylin=linear_fit(y)
    i=0
    while i < len(y):
        y[i] -= ylin[i]
        i += 1
    write_y_to_file(y,"data")

    # Compute the autocorrelation and normalize
    yac=autocorrelate(y)
    yac=normalize_by_first(yac)
    write_y_to_file(yac,"data.ac")

    # The FFT of the autocorrelation
    fac=np.fft.rfft(yac)
    write_abs_y_to_file(fac,"data.fac")

    # Find the (interpolated) dominant mode
    #freq=dominant_freq(fac)
    #period=len(ytypac)/freq
    period=[]
    findperiod=True
    while(findperiod):
        fac=autocorrelate(y)
        yac=normalize_by_first(fac)
        fac=np.fft.rfft(yac)
        p=detect_period(fac,len(yac))

        # Must have 10 cycles in the data.
        if p > len(yac)/10:
            p=1

        if(round):
            p=np.round(p)

        # don't repeat periods
        i=0
        while (i < len(period)):
            if p == period[i]:
                p=1
            i += 1

        period.append(p)
        if p > 1:
            #print p
            #print y[0]
            ytyp=typical_cycle(p,y)
            y=rm_typical_cycle(p,y,ytyp)
            write_y_to_file(ytyp,"data.ytyp"+str(len(period)))
        else:
            findperiod=False
        #if(len(period) > 1):
        #    findperiod=False


        # And stop looping if we get more than 10 (we're probably stuck...)
        if len(period) > 10:
            findperiod=False

    print "Detected period: "+str(period)

    # The data which is not periodic:
    write_y_to_file(y,"data.nac")

    # Write the period length to a file for use with latex.
    f = open('tex/def_period.tex', 'w')
    f.write("\def\periodlength{"+str(period)+"}")
    f.close();
    # Write number of periods to file
    f = open('nperiods', 'w')
    f.write(str(len(period)))
    f.close();

    # Determine the typical cycle:
    #ytyp=typical_cycle(period[0],y)
    #write_y_to_file(ytyp,"data.typ")

    # Determine the part of the signal not represented by the detected
    # cycle:
    #typdiff=typical_cycle_error(period[0],data,ytyp)
    #write_y_to_file(typdiff,"data.dif")    

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
