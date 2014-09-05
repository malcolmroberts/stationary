#!/usr/bin/python -u

import csv # to load the input file.
import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

def stationary_part(list):
    # TODO: actually do some tests here.
    return list

def y_part(list):
    # copy the input list to a buffer for FFTs
    y=[]
    i=0
    while i < len(list):
        y.append(float(list[i][1]))
        i += 1
    return y

def autocorrelate(y):
    ypad=[]
    i=0
    while i < len(y):
        ypad.append(y[i])
        i += 1
    while i < len(y):
        ypad.append(0.0)
        i += 1
    Y=np.fft.rfft(ypad)
    i=0
    while i < len(Y):
        Y[i] *= np.conj(Y[i])
        i += 1
    yac = np.fft.irfft(Y)
    end=np.floor(len(yac)/2)
    yac = yac[0:end]
    norm=yac[0]
    i=0
    while i < len(yac):
        yac[i] /= norm
        i += 1
    return yac

def main(argv):

    usage="./stationary -f <filename>"

    filename=""

    # Load the command-line arguments
    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filename=arg

    # Check that the file exists (and is not a directory)
    if not (os.path.isfile(filename)):
        print "Error: the specified file \""+filename+ " \" does not exist"
        print usage
        sys.exit(2)
    
    # Read the data
    a = []
    csvReader = csv.reader(open(filename, 'rb'), delimiter='\t')
    for row in csvReader:
        a.append(row)
    
    y=y_part(a)
    #print y
    yac=autocorrelate(y)
    print yac

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
