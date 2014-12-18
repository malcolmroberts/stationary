import sys
import numpy as np
import random

# http://www.itl.nist.gov/div898/handbook/eda/section3/eda35d.htm

# NB: this function is not actually used for the Wald-Wolfowotiz test
def count_multiple_runs(y):
    n = len(y)
    if(n == 0):
        return []

    vals = []
    vals.append(y[0])
    last = 0
    runs = []
    runs.append(1)
    
    i = 1
    while(i < n):
        if(not y[i] == vals[last]):
            # see if the previous value is there:
            j = 0
            while(j < len(vals)):
                if(y[i] == vals[j]):
                    break
                j += 1
            if(j < len(vals)):
                runs[j] = runs[j] + 1
            else:
                vals.append(y[i])
                runs.append(1)
            last = j
        i += 1

    return runs, vals

def highlow(y, val):
    highlow = []
    i = 0
    while(i < len(y)):
        if(y[i] > val):
            highlow.append(1)
        if(y[i] < val):
            highlow.append(-1)
        if(y[i] == val):
            highlow.append(2 * random.randint(0, 1) - 1)
        i += 1
    return highlow

def updown(y):
    if(len(y) == 0):
        return []
    if(len(y) == 1):
        return [2 * random.randint(0, 1) - 1]
    updown = []
    i = 1
    while(i < len(y)):
        if(y[i] > y[i-1]):
            updown.append(1)
        if(y[i] < y[i-1]):
            updown.append(-1)
        if(y[i] == y[i-1]):
            updown.append(2 * random.randint(0, 1) - 1)
        i += 1
    return updown

def count_uniques(y):
    counts = []
    vals = []
    i = 0
    while(i < len(y)):
        j = 0
        while(j < len(vals)):
            if(y[i] == vals[j]):
                break
            j += 1
        if(j >= len(vals)):
            vals.append(y[i])
            counts.append(0)
        counts[j] +=  1
        i += 1
    return counts, vals

def count_runs(y):
    if(len(y) == 0):
        return 0
    nruns = 1
    i = 1 
    while(i < len(y)):
        if(y[i] != y[i-1]):
            nruns += 1
        i += 1
    return nruns

def wald_wolfowitz(y):
    counts, vals = count_uniques(y)
    if(len(counts) != 2):
        print "Error : there must be exactly two values in the sequence!"
        sys.exit(1)
    N0 = counts[0]
    N1 = counts[1]
    if(N0 <= 10 or N1 <= 10):
        print "Caution! N0 = " + str(N0) + " and N1 = " +str(N1)\
            + " should both be larger than 10 for statistical validity."
    runs = count_runs(y)
    N = len(y)
    mean = 1 + 2 * N0 * N1 / N
    variance = np.sqrt((mean - 1) * (mean + 1) / (N - 1))
    Z = (runs - mean) / variance
    return Z
