# Do stationarity tests on a list or (time,value) pairs and return the
# stationary tail.
import random
import numpy as np

import scipy.stats # Science!

from utils import *

def stationary_part(list):
    y=y_part(list)

    minlen=100  # min length of sample.

    a=0
    alast=a
    n=len(y)

    # Use bisection method to find start of stationarity
    go=True
    while(go):
        if not 2*np.floor((n-a)/2) == n-a:
            a += 1

        # Not enough points in range: considered not stationary
        if (n - a) < minlen:
            return -1

        h = int(np.floor((n-a)/2))

        if is_stationary(y[n - 2*h:n-h],y[n - h  :n]):
            a = int(np.floor((a+alast)/2))
        else:
            alast=a
            a = int(np.floor((n-a)/2))
        
        # Close enough!
        if(np.abs(a-alast) < 10):
            return a

        # Stop if the entire time series is stationary
        if(a <= 1):
            return a
    
    return a

def is_stationary(y0,y1):
    random.shuffle(y0)
    random.shuffle(y1)

    #T,p= scipy.stats.wilcoxon(y0, y1)

    D,p = scipy.stats.ks_2samp(y0,y1)

    if (p > 0.1):
        # Null hypothesis likely: same dist, so steady 
        return True
    else: 
        return False
