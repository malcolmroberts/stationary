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

# return an array containing the linear fit of the input array of
# values y.
def linear_fit(y):
    A=0.0
    B=0.0
    C=0.0
    D=0.0
    i=0
    while i < len(y):
        A += i
        B += y[i];
        C += i*i;
        D += i*y[i];
        i += 1
    m=(len(y)*D-A*B)/(len(y)*C-A*A);
    c=(B-m*A)/len(y);
    #print "m="+str(m)
    #print "c="+str(c)
    ylin=[]
    i=0
    while i < len(y):
        ylin.append(m*i +c)
        i += 1
    return ylin

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

# Return the absolute value of the complex value z: 
def abs(z):
    return np.sqrt(z.real*z.real + z.imag*z.imag)

# Return the max of the quadratic spline going through three equally
# spaced points with y-values y0, y1, and y2.
def paramax(y0, y1, y2):
    return 0.5*(y0 - y2)/(y0 - 2*y1 +y2)

# Return the index of the largest mode in a Fourier series Y, the
# interpolate (using a quadratic approximation around the max) to find
# the actual max.
def dominant_mode(Y):
    max=0.0
    i=0
    imax=0
    while i < len(Y):
        amp=abs(Y[i])
        if amp > max:
            max=amp
            imax=i
        i += 1
    return imax+paramax(abs(Y[imax-1]),abs(Y[imax]),abs(Y[imax+1]))

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

def typical_period(period, y):
    # The length of the the output array is the floor of the period length:
    n=int(np.floor(period))

    # Number of periods in y:
    nperiod=int(np.floor(len(y)/period))
    print "There were "+str(nperiod)+" period(s) found in the data."
    
    ytyp=[]
    i=0
    while i < n:
        ytyp.append(0.0)
        j=0
        while j < nperiod:
            jbase=int(np.floor(period*j))
            ytyp[i] += y[jbase+i]
            j += 1
        ytyp[i] /= nperiod
        i += 1
    return ytyp

# Main program
def main(argv):
    usage="./stationary -f <filename>"

    filename=""

    # Load the command-line arguments:
    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filename=arg

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

    # Output the data with linear regression removed:
    datawriter = csv.writer(open("data", 'wb'), delimiter='\t')
    i=0
    while i < len(y):
        datawriter.writerow([i,y[i]])
        i += 1

    # Compute the autocorrelation and normalize:
    yac=autocorrelate(y)
    yac=normalize_by_first(yac)
    
    # Output the autocorrelation:
    datawriter = csv.writer(open("data.ac", 'wb'), delimiter='\t')
    i=0
    while i < len(yac):
        datawriter.writerow([i,yac[i]])
        i += 1

    # The FFT of the autocorrelation:
    fac=np.fft.rfft(yac)

    # Output the FFT of the audotorrelation:
    datawriter = csv.writer(open("data.fac", 'wb'), delimiter='\t')
    i=0
    while i < len(fac):
        datawriter.writerow([i,abs(fac[i])])
        i += 1

    # Find the (interpolated) dominant mode:
    freq=dominant_mode(fac)
    period=len(yac)/freq
    print "Detected period: "+str(period)

    # Determine the typical cycle:
    ytyp=typical_period(period,y)

    # Output the typical period:
    datawriter = csv.writer(open("data.typ", 'wb'), delimiter='\t')
    i=0
    while i < len(ytyp):
        datawriter.writerow([i,ytyp[i]])
        i += 1

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
