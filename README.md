# calc_hdx

calc_hdx is a Python 3 package used to analyze molecular dynamics trajectories and calculate predicted HDX-MS protection factors and deuterated fractions. The predictions are optionally compared with experimental data, and results plotted in a series of output pdf files. The code currently (v0.1) contains the following files:

- calc_hdx.py : Main executable, usable from the command line with provided arguments, or by importing the calc_hdx package

- Functions.py : Transferable functions for trajectory reading/manipulating topologies and trajectory analysis

- DfPred.py : Parent class for analysis methods. Functions for calculating intrinsic rates from a topology, and for calculating deuterated fractions

- Methods.py : Analysis method classes for calculating protection factors with various models. Currently contains 'Radou' model (identical to Best & Vendruscolo, Structure, 2006, 14 (1), 97-106), and 'Persson-Halle' model (Persson & Halle, PNAS, 2015, 112 (33), 10383-10388).

- Analysis.py : Classes for calculating deuterated fractions and outputting results in tabular and graphical form

## Dependencies & usage instructions

In addition to the standard Python 3.7 libraries it has been developed with the following dependencies:

MDtraj 1.9.3
Numpy 1.16.4
Scipy 1.3.0
Matplotlib 3.1.0

No testing has been performed with earlier or later dependency versions.

For a description of calc_hdx command line options, run `calc_hdx.py -h`

For an example usage and dataset, see DOI: 10.5281/zenodo.3385169

