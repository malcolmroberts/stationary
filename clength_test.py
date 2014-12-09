#!/usr/bin/python -u

from cycle import *
from utils import *
from stats import *
import random

def ar1(n, theta):
    y = []
    y0 = 0
    for i in range(n):
        y.append(theta * y0 + random.uniform(-1,1))
        y0 = y[i]
    return y

n = 10000
N = 100
theta = 0.0

meanlength = 0.0

i = 0
while i < N:
    y = ar1(n, theta)
    meanlength += correlation_length(y)
    i += 1

print "mean correlation length:"
print(meanlength / N)
