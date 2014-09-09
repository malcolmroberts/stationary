#!/usr/bin/python -u

import csv # to load the input file.
import numpy # numerical methods such as FFTs.
import sys # to check for file existence, etc.
import getopt # to pass command-line arguments
import os.path #for file-existence checking

def main(argv):
    usage="./change_format.py <input filename> <output filename>"
    
    if len(sys.argv) < 2:
        print "Please specify an input."
        print usage
        exit(1)

    infile=argv[0]

    outfile=infile+".out"

    if len(sys.argv) == 3 :
        outfile=argv[1]

    # Check that the file exists (and is not a directory)
    if not (os.path.isfile(infile)):
        print "Error: the specified file \""+infile+ " \" does not exist"
        print usage
        sys.exit(2)
    
    # Read the data
    a = []
    csvReader = csv.reader(open(infile, 'rb'), delimiter='\t')
    for row in csvReader:
        a.append(row)

    # reformat the data and save in b.
    b = []
    i=0
    while (i < len(a)):
        b.append([i, a[i][0]])
        i += 1

    # write the data out.
    with open(outfile, 'w') as fp:
        B = csv.writer(fp, delimiter='\t')
        B.writerows(b)

# The main program is called from here
if __name__ == "__main__":
    main(sys.argv[1:])
