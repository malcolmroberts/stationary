# Do stationarity tests on a list or (time,value) pairs and return the
# stationary tail.

import random
import numpy as np
import sys

# From the current repo:
from cycle import *
from wald_wolfowitz import *

import scipy.stats # Science!

from utils import *

def bin_sequence(y, binlen):
    nbins = len(y) / binlen
    ybins = []
    i = 0
    while(i < nbins):
        ybins.append(0)
        j = 0
        while(j < binlen):
            ybins[i] += y[i * binlen + j]
            j += 1
        ybins[i] /= binlen
        i += 1
    return ybins

# Return the correlation length of the sequence of floats y.
def correlation_length(y):
    yac = autocorrelate(y)
    yac = normalize_by_first(yac)
    n = len(yac)

    # 95% confidence interval for autocorrelation
    ac95 = 1.96 / np.sqrt(n)

    wittetal = False

    firstcross = False

    acarea = False

    acabsarea = False

    totacarea = False

    totacabsarea = True

    if(wittetal):
        i = n
        while i >= 0:
            i -= 1
            #ac95 = (-1 + 1.96 * np.sqrt(n - i - 1)) / (n-i)
            if np.abs(yac[i]) > ac95:
                return i
        return 0

    if(firstcross):
        i = 1
        while i < n:
            ac95 = (-1 + 1.96 * np.sqrt(n - i - 1)) / (n-i)
            if yac[i] < ac95:
                return i - 1
            i += 1

    if(acarea):
        area = 0
        i = 1
        while(i < n and np.abs(yac[i]) < ac95):
            area += yac[i]
            i += 1
        return area / 1.0

    if(acabsarea):
        area = 0
        i = 1
        while(i < n and np.abs(yac[i]) < ac95):
            area += abs(yac[i])
            i += 1
        return area / 1.0

    if(totacarea):
        area = 0
        i = 1
        while(i < n):
            area += yac[i]
            i += 1
        return area / 1.0

    if(totacabsarea):
        area = 0
        i = 1
        while(i < n):
            area += abs(yac[i])
            i += 1
        return area / 1.0

    return n

def stationary_part(list, stest, p, rmcycles, roundperiod):
    y = y_part(list)

    minlen = 500  # min length of sample.
    n = len(y)

    a = 0
    low = 0
    high = n

    # Use bisection method to find start of stationarity
    go = True
    while(go):
        # Not enough points in range: considered not stationary
        if (n - a) < minlen:
            return -1

        alast = a
        ytest = y[a:n]

        ytestpow = power(ytest)
        if(rmcycles == True):
            cycles, ytest = find_multiple_periods(ytest, roundperiod)
        # If the non-periodic part of the signal is negligible,
        # consider the signal stationary
        if(power(ytest) / ytestpow < 1e-12):
            return a
        
        good = is_stationary(ytest, stest, p, minlen)

        if good:
            high = a
            a = int(np.floor((a + low) / 2))
        else:
            low = a
            a = int(np.floor((high + a) / 2))

        # Close enough!
        if(np.abs(a - alast) < 2):
            return a

        # Stop if the entire time series is stationary
        if(a <= 1):
            return a
    
    return a

def is_stationary(y, stest, p_crit, minlen):
    n = len(y)
    p = -1

    if stest == "wsr":
        h = int(np.floor(n / 2))
        y0 = y[0:h]
        y1 = y[n-h:n]
        
        random.shuffle(y0)
        random.shuffle(y1)

        repeats, nr = scipy.stats.find_repeats(y)
        ni = n
        i = 0
        while i < len(nr):
            ni -= nr[i]
            i += 1
        if ni < minlen:
            p = 1 # we call this stationary
        else :
            T , p = scipy.stats.wilcoxon(y0, y1)

    if stest == "ks":
        h = int(np.floor(n / 2))
        y0 = y[0:h]
        y1 = y[n-h:n]
        D, p = scipy.stats.ks_2samp(y0, y1)

    if stest == "highlowruns" or stest == "updownruns":
        clen = correlation_length(y)
        ybins = bin_sequence(y, int(round(clen)) + 1)
        if(len(ybins) < 2):
           return False

        x = []
        if stest == "highlowruns":
            mean = float(sum(ybins)) / len(ybins)
            x = highlow(ybins, mean) # is the mean a good fit for the data?

        if stest == "updownruns":
            x = updown(ybins) # is the mean a good fit for the data?

        counts, vals = count_uniques(x)
        if(len(vals) == 1):
            return True # a constant sequence is stationary
        if(len(vals) != 2):
            print "Error: runs test is considering " + str(len(vals)) \
                + " different values: cannore be more than two!"
            print vals
            exit(1)
        if(max(counts[0], counts[1]) < 10):
            return False

        Z = wald_wolfowitz(x)
        if(np.abs(Z) > 1.96):
            return False
        else:
            return True

    # requires 0.14 of scipy
    #A2, critical, p = scipy.stats.anderson_ksamp([y0,y1])
    
    if(p == -1):
        print "Error: invalid choice of statistical test "+stest
        exit(1)

    if (p > p_crit):
        # Null hypothesis likely: same dist, so steady 
        return True
    else: 
        return False

def bindata(y, N):
    z = []
    n = int(np.floor(len(y) / N)) # number of bins
    i = 0    
    while(i < n):
        z.append(0.0)
        j = 0
        while(j < N):
            z[i] += y[N * i + j]
            j += 1
        z[i] /= N
        i += 1
    return z


def variance(z):
    zbar = 0
    i = 0
    while(i < len(z)):
        zbar += z[i]
        i += 1
    zbar /= len(z)
    var = 0
    i = 0
    while(i < len(z)):
        d = zbar - z[i]
        var += d * d
        i += 1
    var /= len(z)
    return var

def mser5(y):
    N = 5 # size of bins
    z = bindata(y, N)
    dstar = -1
    vmin = sys.float_info.max
    nw = len(z)
    d = 0
    while(d < len(z)):
        val = variance(z[d:nw]) / (nw - d)
        if(val < vmin):
            vmin = val
            dstar = d
        d += 1
    print vmin
    return dstar * N


def mae(z):
    zbar = 0
    i = 0
    while(i < len(z)):
        zbar += z[i]
        i += 1
    zbar /= len(z)
    var = 0
    i = 0
    while(i < len(z)):
        d = np.abs(zbar - z[i])
        var += d
        i += 1
    var /= len(z)
    return var

def maer5(y):
    N = 5 # size of bins
    z = bindata(y, N)
    dstar = -1
    vmin = sys.float_info.max
    nw = len(z)
    d = 0
    while(d < len(z)):
        val = mae(z[d:nw]) / (nw - d)
        if(val < vmin):
            vmin = val
            dstar = d
        d += 1
    print vmin
    return dstar * N

