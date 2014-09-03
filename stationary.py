#!/usr/bin/python -u

import csv # to load the input file.
import numpy # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

def main(argv):

    usage="FIXME"

    filename=""

    # Load the command-line arguments
    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            filename=arg

    # Check that the file exists (and is not a directory)
    if not (os.path.isfile(filename)):
        print "Error: the specified file \""+filename+ " \" does not exist"
        print usage
        sys.exit(2)
    
    # Read the data
    a = []
    csvReader = csv.reader(open(filename, 'rb'), delimiter='\t')
    for row in csvReader:
        a.append(row)
    
    

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
