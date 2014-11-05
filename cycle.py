import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.

# Return an array containing the linear fit of the input array of
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
    maxit=64 # Maximum number of iterations
    maxerr=1e-3 # Tolerance for the difference between the detected
                # mode and the nearest mode on the grid.
    
    # Set test-max to a very large value initially
    err=sys.float_info.max

    # Initialize nfit
    nfit=n

    # nfit0 is used to see if nfit is stable from iteration to
    # iteration.
    nfit0=0

    cont=True
    i=0
    while cont:
        # Mode (on the grid) with max frequency:
        freq_i = dominant_freq(Y[0:int(np.floor(nfit/2))])

        if freq_i == 0:
            return 1

        # Interpolate around the peak to see if there's a better max:
        ishift=paramax(abs(Y[freq_i-1]),abs(Y[freq_i]),abs(Y[freq_i+1]))
        # Interpolated frequency:
	freq=freq_i + ishift

        # Length of period from interpolated frequency:
        period=n/freq
        
        # Number of interpolated periods that fit in data:
        nperiod=int(np.floor(n/period))
        
        # Max data size with an integral number of periods:
        nfit=int(np.floor(nperiod*period))
        
        # Difference between frequency and nearest integral mode:
        err=np.abs(freq-np.round(freq))

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
    nperiod=int(np.ceil(len(y)/period))
    i=0
    while i < len(ytyp):
        j=0
        while j < nperiod:
            jbase=int(np.floor(period*j))
            pos=jbase+i
            if(pos < len(y)):
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


def significant_cycle(yac):
    izero=0

    i=0
    while(izero == 0 and i < len(yac)):
        if(yac[i] < 0):
            izero=i
        i += 1
        
    n=len(yac)

    i=izero
    while(i < n):
        #ac95 = 1.96 / np.sqrt(n)
        ac95 = (-1 + 1.96 * np.sqrt(n - i - 1))/(n-i)
        if(np.abs(yac[i]) > ac95):
            return True
        i += 1

    return (False)


# Find all periods in the input data y.
# Input: sequence y, bool round (if we're rounding to integral periods)
# Output: cycles, which contains [period length , sequence of typical period]
# with the last element being [1, sequence with non-periodic part]
def find_multiple_periods(y,round):
    cycles=[]
    findperiod=True
    while(findperiod):
    
        yac=autocorrelate(y)
        yac=normalize_by_first(yac)

        if not significant_cycle(yac):
           return cycles

        fac=np.fft.rfft(yac)
        p=detect_period(fac,len(yac))
    
        # Must have 10 cycles in the data.
        if p > len(yac)/10:
            p=1
        
        if p > 1:
            if(round):
                p=np.round(p)

        # don't repeat periods
        i=0
        while (i < len(cycles)):
            if p == cycles[i][0]:
                p=1
            i += 1

        if p > 1:
            ytyp=typical_cycle(p,y)
            cycles.append([p,ytyp])
            y=rm_typical_cycle(p,y,ytyp)
        else:
            findperiod=False

        # And stop looping if we get more than 10 (we're probably stuck...)
        if len(cycles) > 20:
            findperiod=False
            
    cycles.append([1,y])
    return cycles, y
