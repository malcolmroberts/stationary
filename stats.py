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

    # Play high-low
    go=True
    while(go):
        if not 2*np.floor((n-a)/2) == n-a:
            a += 1

        if n - a < minlen:
            return -1

        h = int(np.floor((n-a)/2))
        y0=y[n - 2*h:n-h]
        y1=y[n - h  :n]
        s=is_stationary(y0,y1)

        if s:
            a = int(np.floor((a+alast)/2))
        else:
            alast=a
            a = int(np.floor((n-a)/2))
        
        if(np.abs(a-alast) < 10):
            return a

        if(a <= 1):
            return a
    
    return a

def is_stationary(y0,y1):
    random.shuffle(y0)
    random.shuffle(y1)
    T,p= scipy.stats.wilcoxon(y0, y1)
    #print "T="+str(T)
    #print "p="+str(p)
    if (p > 0.1):
        # Null hypothesis likely: same dist, so steady 
        return True
    else: 
        return False
