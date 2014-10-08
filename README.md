stationary
==========

A Python script to determine when a given time series has reached a
statistically stationary state and if the final signal is periodic.

Syntax:

    ./stationary.py -f <FILE> 

where the input file FILE has the index (as an integer) followed by a
the value of the signal (as a float) separated by a tab.  If the data
is not in this format, it can be converted to this format by running

    ./change_format.py <IN> <OUT>

with the <out> argument optionally specifying the output filename (the
default output filename is <IN>.out).

Use:

stationary.py currently produces a variety of files:

* data, which is the input with the linear regression removed,

* data.ac, which is the autocorrelation of data,

* data.fac, which is the magnituded of the modes of the DFT of the
  autocorrelation,

* data.typ<#>, which are the typical cycles detected (<#> an ineger >= 0).

* data.dif, which is the difference between the signal (with the
  linear regression already removed) and the signal as represented by
  the typical cycle.  The program also outputs (to the terminal) the
  RMS difference between the signal and its represenatation by the
  typical cycle.


Viewing the output:

The output can be used with the included plot.asy file, which is
opened with Asymtote (asymptote.sf.net).  Running plot.asy with the
command

  asy plot.asy

The script then asks for the filenames, which is a comma-separated
list of filenames (eg "data,data.typ", without the quotes), and then
asks for a choice of scales and axis labels.


Scripts:
The file
  run.sh <filename> <optional bool to specify rounding>
take an input file <filename> and convert it into a sequence of
(position,value) pairs, and then run stationary.py

The file
  run2.sh <filename> <optional bool to specify rounding>
calls run.sh and then creates a PDF of the output, using asymptote and
tex for visualization and typesetting.
