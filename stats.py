# Do stationarity tests on a list or (time,value) pairs and return the
# stationary tail.
def stationary_part(list):
    stable=list # TODO: actually do some tests here.
    #stable=list[0:200]
    return stable

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
