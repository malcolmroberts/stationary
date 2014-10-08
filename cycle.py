import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.

# Return the autocorrelation of the array of reals y.
def autocorrelate(y):
    ypad=[]
    i=0
    while i < len(y):
        ypad.append(y[i])
        i += 1
    i=0
    while i < len(y):
        ypad.append(0.0)
        i += 1
    norm=1.0/len(y)
    Y=np.fft.rfft(ypad)
    i=0
    while i < len(Y):
        Y[i] *= np.conj(Y[i])
        i += 1
    yac = np.fft.irfft(Y)
    yac = yac[0:len(y)]
    # Normalize by the length of the FFT to get the bare
    # autocorrelation:
    i=0
    while i < len(yac):
        yac[i] *= norm;
        i += 1
    return yac

# Normalize an array by the value of the first element.
def normalize_by_first(y):
    norm=y[0]
    if(np.abs(norm) > 0.0):
        i=0
        while i < len(y):
            y[i] /= norm
            i += 1
    return y

# Return the max of the quadratic spline going through three equally
# spaced points with y-values y0, y1, and y2.
def paramax(y0, y1, y2):
    if(np.abs(y0) > 0 and np.abs(y1) > 0 and np.abs(y2) > 0): 
        return 0.5*(y0 - y2)/(y0 - 2*y1 +y2)
    return 0

# Return the index of the largest mode in a Fourier series Y, the
# interpolate (using a quadratic approximation around the max) to find
# the actual frequency corresponding to the maximum.
# Y is the FFT of the signal or its autocorrelation.
def dominant_freq(Y):
    max=0.0
    i=0
    imax=0

    while i < len(Y):
        amp=abs(Y[i])
        if amp > max:
            max=amp
            imax=i
        i += 1

    return imax



# Return the dominant frequency using an iterative algorithm:
# Y is the FFT of the signal or its autocorrelation.
# n is the length of the original data.
def detect_period(Y,n):
    maxit=32 # Maximum number of iterations
    maxerr=1e-3 # Tolerance for the difference between the detected
                # mode and the nearest mode on the grid.
    
    # Set test-max to a very large value initially
    err=sys.float_info.max
    # Initialize n
    nfit=n
    # nfit0 is used to see if nfit is stable from iteration to
    # iteration.
    nfit0=0

    cont=True
    i=0
    while cont:
        # Mode (on the grid) with max frequency:
        freq_i=dominant_freq(Y[0:int(np.floor(nfit/2))])
        #print freq_i
        #print len(Y)

        # Interpolate around the peak to see if there's a better max:
        ishift=paramax(abs(Y[freq_i-1]),abs(Y[freq_i]),abs(Y[freq_i+1]))
        # Interpolated frequency:
	freq=freq_i + ishift

        # Amplitude of dominant mode of the autocorrelation:
        amp=abs(Y[freq_i])
        #print "amplitude="+str(amp)
        # Confidence interval:
        nf=1.96/np.sqrt(n);
        #print "Confidence interval: "+str(nf)

        # check whether the amplitude is above the confidence interval:
        if amp < nf:
            print "amp too small"
            return 1

        # Length of period from interpolated frequency:
        period=n/freq
        # Number of interpolated periods that fit in data:
        nperiod=int(np.floor(n/period))
        # Max data size with an integral number of periods:
        #nfit=int(np.floor(nperiod*period)/2)
        nfit=int(np.floor(nperiod*period))
        # Difference between frequency and nearest integral mode:
        err=np.abs(freq-np.round(freq))
        
        # Diagnostic output:
        #print
        #print "nfit="+str(nfit)
        #print "nfit0="+str(nfit0)
        #print "freq="+str(freq)
        #print "period="+str(period)
        #print "nperiod="+str(nperiod)
        #print "err="+str(err)
       
        # Stop if the number of cycles that fits changes by less than
        # 1 data point per iteration.
        if (np.fabs(nfit-nfit0) < 1):
            #print ("Finisehd after "+str(i)+" iterations: nfit="
            #       + str(nfit)+" is stable.")
            cont=False
        nfit0=nfit
        
        # Stop if the frequency falls close to a Fourier mode.
        if (err < maxerr):
            #print "Done: error is small: "+str(err)+"<"+str(maxerr)
            cont=False

        if (i > maxit):
            #print "Did "+str(i)+" iterations: stopping."
            
            #print "nfit="+str(nfit)
            #print "nfit0="+str(nfit0)
            #print "freq="+str(freq)
            #print "period="+str(period)
            #print "nperiod="+str(nperiod)
            #print "err="+str(err)
            cont=False

        i += 1

    return period


# Return an array containing the typical cycle of the signal.
# period is the period of the cycle.
# y is the original signal
def typical_cycle(period, y):
    # The length of the the output array is the floor of the period length:
    n=int(np.floor(period))

    # Number of periods in y:
    nperiod=int(np.floor(len(y)/period))
    #print "There were "+str(nperiod)+" period(s) found in the data."
    
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
        
    error=0.0
    i=0
    while i < n:
        j=0
        while j < nperiod:
            jbase=int(np.floor(period*j))
            norm=np.fmax(np.fabs(ytyp[i]),np.fabs(y[jbase+i]))+1e-7
            diff=(ytyp[i] - y[jbase+i])/norm
            error +=  diff*diff
            j += 1
        i += 1
    #print "RMS difference: "+str(np.sqrt(error/len(y)))

    return ytyp

# FIXME: add documentation
def rm_typical_cycle(period,y,ytyp):
    newy=[]
    i=0
    while i < len(y):
        newy.append(y[i])
        i += 1
    nperiod=int(np.floor(len(y)/period))
    i=0
    while i < len(ytyp):
        j=0
        while j < nperiod:
            jbase=int(np.floor(period*j))
            pos=jbase+i
            diff=ytyp[i] - float(y[pos])
            newy[pos] = diff
            j += 1
        i += 1
    return newy



# Return and array containing the L2 distance between the typical
# cycle and the actual data.
# Input:
# period is the period (not necessarily an integer)
# data is the original data, in (time,value) pairs
# Output:
# ytyp is an array containing the typical signal. 
def typical_cycle_error(period,data,ytyp):
    typdiff=[]
    i=0
    while i < len(data):
        typdiff.append([data[i][0],0])
        i += 1
    nperiod=int(np.floor(len(data)/period))
    i=0
    while i < len(ytyp):
        j=0
        while j < nperiod:
            jbase=int(np.floor(period*j))
            pos=jbase+i
            diff=ytyp[i] - float(data[pos][1])
            typdiff[pos][1] = diff
            j += 1
        i += 1
    return typdiff
