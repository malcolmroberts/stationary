#!/usr/bin/python -u

import csv # to load the input file.
import numpy as np # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

# Load files in this repository
from utils import *
from stats import *

# Main program
def main(argv):
    usage = "Usage: FIXME\n"
    filename = ""

    # Load the command-line arguments
    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filename = arg

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

    y = y_part(data)
    
    print maer5(y)


# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
