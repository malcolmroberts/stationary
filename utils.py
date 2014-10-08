import numpy as np # numerical methods such as FFTs.

import csv # to load the input file.

# Return the list of values from a list of (time,value) pairs.
def y_part(list):
    y=[]
    i=0
    while i < len(list):
        y.append(float(list[i][1]))
        i += 1
    return y

# Return the square of the absolute value of the complex value z:
def abs2(z):
    return z.real*z.real + z.imag*z.imag

# Return the absolute value of the complex value z: 
def abs(z):
    return np.sqrt(abs2(z))

# write the y-values to a an output file in (t,val) pairs.
def write_y_to_file(y,filename):
    f=open(filename, 'wb')
    datawriter = csv.writer(f, delimiter='\t')
    i=0
    while i < len(y):
        datawriter.writerow([i,y[i]])
        i += 1
    f.close()

# write the y-values to a an output file in (t,val) pairs.
def write_abs_y_to_file(Y,filename):
    f=open(filename, 'wb')
    datawriter = csv.writer(f, delimiter='\t')
    i=0
    while i < len(Y):
        datawriter.writerow([i,np.abs(Y[i])])
        i += 1
    f.close()
