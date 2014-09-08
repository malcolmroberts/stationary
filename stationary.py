#!/usr/bin/python -u

import csv # to load the input file.
import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

# Do stationarity tests on a list or (t,value) pairs and return the
# stationary tail.
def stationary_part(list):
    stable=list # TODO: actually do some tests here.
    return stable

# Return the list of values from a list of (t,value) pairs.
def y_part(list):
    y=[]
    i=0
    while i < len(list):
        y.append(float(list[i][1]))
        i += 1
    return y

# Return the autocorrelation of the array of reals y.
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
    return yac

# Normalize an array by the value of the first element.
def normalize_by_first(y):
    norm=y[0]
    i=0
    while i < len(y):
        y[i] /= norm
        i += 1
    return y

# Return the index of the largest mode in a Fourier series Y.
def dominant_mode(Y):
    max=0.0
    i=0
    imax=0
    while i < len(Y):
        amp=Y[i].real*Y[i].real + Y[i].imag*Y[i].imag
        if amp > max:
            max=amp
            imax=i
        i += 1
    return imax

# Return the max of the quadratic spline going through three equally
# spaced points with y-values y0, y1, and y2.
def paramax(y0, y1, y2):
    0.5*(y0 - y2)/(y0 - 2*y1 +y2)

# Return the value of the quadratic spline at position x between
# equally spaced points with input values y0, y1, and y2.
def spline(x, y0, y1, y2):
    x0=-1.0
    x1=0.0
    x2=1.0
    L0=1.0
    L0 *= (x-x1)/(x0-x1)
    L0 *= (x-x2)/(x0-x2)

    L1=1.0
    L1 *= (x-x0)/(x1-x0)
    L1 *= (x-x2)/(x1-x2)
    
    L2=1.0
    L1 *= (x-x0)/(x2-x0)
    L0 *= (x-x1)/(x2-x1)
  
    return y0*L0 + y1*L1 + y2*L2

# Main program
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
    # TODO: remove the linear fit (or just the mean?)
    #print y
    yac=autocorrelate(y)
    yac=normalize_by_first(yac)
    #print yac

    fac=np.fft.rfft(yac)

    #print fac
    print dominant_mode(fac)


# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
