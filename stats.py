# Do stationarity tests on a list or (time,value) pairs and return the
# stationary tail.
import random
import numpy as np
from cycle import *

import scipy.stats # Science!

from utils import *

def stationary_part(list, stest, p, rmcycles, roundperiod):
    y = y_part(list)

    minlen = 100  # min length of sample.
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
         
        # if the non-periodic part of the signal is negligible,
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

    #cycles, y=find_multiple_periods(y,round)

    h = int(np.floor(n / 2))
    y0 = y[0:h]
    y1 = y[n-h:n]

    random.shuffle(y0)
    random.shuffle(y1)

    p = -1

    if stest == "wsr":
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
        D, p = scipy.stats.ks_2samp(y0, y1)

    # requires 0.14 of scipy
    #A2, critical, p = scipy.stats.anderson_ksamp([y0,y1])
    
    if(p == -1):
        print "Error: invalid choice of statistical test "+stest
        exit(1)


    #print p

    if (p > p_crit):
        # Null hypothesis likely: same dist, so steady 
        return True
    else: 
        return False
