#!/usr/bin/python -u

import csv # to load the input file.
import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

from cycle import *
from utils import *
from stats import *
    
def process_stationary_signal(start,data,round):
    # Put the y-values from the input data into y:
    y = y_part(data)

    # Remove the mean signal:
    ymean = 0
    i = 0
    while i < len(y):
        ymean += y[i]
        i += 1
    i = 0
    while i < len(y):
        y[i] -= ymean / len(y)
        i += 1

    # # Remove the linear fit:
    # ylin=linear_fit(y)
    # i=0
    # while i < len(y):
    #     y[i] -= ylin[i]
    #     i += 1

    # Compute the autocorrelation and normalize
    yac = autocorrelate(y)
    yac = normalize_by_first(yac)
    #write_y_to_file(yac, "data.ac")

    # The FFT of the autocorrelation
    fac=np.fft.rfft(yac)
    write_abs_y_to_file(fac, "data.fac")

    # Find all periodic components, store in the cycles array.
    # Last element of array contains remainder.
    cycles, yleft = find_multiple_periods(y, round)

    return yac, fac, cycles, yleft


# Main program
def main(argv):
    usage = "Usage:\n"\
        "./stationary\n"\
        "\t-f <filename> input filename\n"\
        "\t-r <0 or 1 (default)>  round period to nearest ineteger?\n"\
        "\t-s <ks (default), wsr>\ statistical test\n"\
        "\t-p <real, defaul=0.1>\ specify p-value for stationarity\n"\
        "\t-c <0 or 1 (default)>  Remove period before testing for stationairty\n"\
        "\t-h display help message\n"

    # Filename for the input data (consisting of tab-separated
    # (time,value) pairs.
    filename=""

    preremove_cycles = True
    round = True
    #round=False
    p = 0.1
    stest="ks"

    # Load the command-line arguments
    try:
        opts, args = getopt.getopt(argv,"c:f:r:p:s:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-p"):
            p = float(arg)
        if opt in ("-s"):
            stest=arg
        if opt in ("-f"):
            filename = arg
        if opt in ("-r"):
            round = (arg == "True" or arg == "true" or arg == "1")
        if opt in ("-c"):
            preremove_cycles = (arg == "True" or arg == "true" or arg == "1")

    print "Using "+stest+ " statistical test with p="+str(p)

    if round:
        print "Period lengths are rounded."
    else:
        print "Period lengths are not rounded."

    if preremove_cycles:
        print "Cyclic component removed before stationarity test."
    else:
        print "Cyclic component included in stationarity test."

    print

    # Check that the file exists (and is not a directory):
    if not (os.path.isfile(filename)):
        print "Error: the specified file \"" + filename + " \" does not exist"
        print usage
        sys.exit(2)
    else:
        print "Analyzing " + filename

    # Read the data
    data = []
    csvReader = csv.reader(open(filename, 'rb'), delimiter = '\t')
    for row in csvReader:
        data.append(row)
    print str(len(data)) + " data points found"

    # Write the input to file.
    write_tv_seq_to_file(data,"data.in")

    ypower = power(y_part(data))
    print "power of input: " + str(ypower)

    corrlen = correlation_length(y_part(data))
    print "Correlation length of input is " + str(corrlen)

    # Consider only the stationary part of the data:
    start = stationary_part(data, stest, p, preremove_cycles, round)

    # Write the period length to a file for use with latex.
    f = open('tex/def_start.tex', 'w')
    f.write("\def\startval{" + str(start) + "}")
    f.close();

    # Output the start of stationarity for the bash script
    f = open('startval', 'w')
    f.write(str(start))
    f.close();

    # Output for creating the summary pdf
    f = open("output", 'wb')
    f.write(str(ypower))
    f.write("\t")
    f.write(str(start))
    f.write("\t")
    f.close()

    if(start < 0):
        print "Signal is non-stationary."
    else:
        print "Stationarity starts at " + str(start) +  " of " \
            + str(len(data)) + " points (" + str((100.0*start) / len(data))\
            +" %)"
        data = data[start:len(data)]
        write_tv_seq_to_file(data,"data.stat")
        yac, fac, cycles, yleft = process_stationary_signal(start,data,round)
    
        corrlen = correlation_length(yleft)
        print "Correlation length of remainder is " + str(corrlen)

        write_y_to_file(yac, "data.ac")
        write_abs_y_to_file(fac, "data.fac")

        i = 0
        while i < len(cycles)-1:
            print "Period length " + str(cycles[i][0]) \
                + " has power " + str(power(cycles[i][1]))
            i += 1
        print "Non-periodic part of signal has power "\
            + str(power(cycles[i][1]))

        # Compute the autocorrelation and normalize
        write_y_to_file(yleft,"data.np",start)

        #print "Detected period: "+str(period)
        periods = []
        i = 0
        while i < len(cycles)-1:
            #print "period="+str(cycles[i][0])
            periods.append(cycles[i][0])
            write_y_to_file(cycles[i][1], "data.ytyp" + str(i))
            i += 1

        # Output the number of periods for the bash script
        f = open('nperiods', 'w')
        f.write(str(len(periods)))
        f.close();

        # Output the number of periods for the tex file
        f = open('tex/def_nperiods.tex', 'w')
        f.write("\def\\nperiods{" + str(len(periods)) + "}")
        f.close();
    
        if len(periods) > 0:
            # Write the period length to a file for use with latex.
            f = open('tex/def_period.tex', 'w')
            f.write("\def\periodlength{" + str(periods) + "}")
            f.close();
    
            periodpower = []
            i = 0
            while i < len(periods):
                periodpower.append(power(cycles[i][1]))
                i += 1
            f = open('tex/def_period_power.tex', 'w')
            f.write("\def\periodpower{" + str(periodpower) + "}")
            f.close();
        
        # Write number of periods to file
        f = open('nperiods', 'w')
        f.write(str(len(periods)))
        f.close();

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
