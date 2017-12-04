#!/usr/bin/env python

# Wrapper script for HDX analyses
# Author: Richard Bradshaw, richard.bradshaw@nih.gov
#
# For help/usage instructions: calc_hdx.py -h
#
#
# Python 3 compatibilities
from __future__ import print_function
from __future__ import division
# Dependencies
import mdtraj as md
import sys, ast
import argparse
# 
import Functions, Methods


### Argparser ###
def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--traj",help="Trajectory/ies for analysis",nargs='+',type=str, required=True)
    parser.add_argument("-p","--parm",help="Topology file to be used for analysis",type=str, required=True)
    parser.add_argument("-s","--stride",help="Stride at which to read the trajectory. Default = 1 (every frame)", nargs=1, type=int, default=1)
    parser.add_argument("-c","--chunks",help="If set, trajectory will be read in chunks of this size (lowers memory requirements for large trajectories). Default = 1000", nargs='?', type=int, const=1000)
    parser.add_argument("-m","--method",help="Method for analysis. Currently 'Radou' is the only option", choices=['Radou'], default='Radou', required=True)
    parser.add_argument("-dt","--times",help="Times for analysis, in minutes. Defaults to [ 0.167, 1.0, 10.0, 120.0 ]", nargs='+', default=[0.167, 1.0, 10.0, 120.0])
    parser.add_argument("-log","--logfile",help="Name of logfile for printout of run info. Defaults to 'HDX_analysis.log'", type=str, default='HDX_analysis.log')
    parser.add_argument("-seg","--segfile",help="Name of file with segment definitions for analysis. Segments should be defined one per line, with starting/finishing residues whitespace separated. Defaults to 'cropped_seg.list'",type=str, default='cropped_seg.list')
    parser.add_argument("-out","--outprefix",help="Prefix for prediction output files",type=str, default='')
    parser.add_argument("-opt","--method_options",help="Additional method options. Should be provided as a single string in Python dictionary format, e.g.:  '{ 'hbond_method' : 'contacts', 'cut_Nc' : 0.70, 'save_contacts' : True }' (Note the string must be enclosed in quotes)",type=str)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    if args.method_options is not None:
        try:
            optdict = ast.literal_eval(args.method_options)
            if isinstance(optdict, dict):
                args.method_options = optdict
            else:
                raise Functions.HDX_Error("Your options flag isn't a dictionary. Dictionary format with key/value pairs is required")
        except ValueError:
            raise Functions.HDX_Error("There's something wrong with the syntax of your options flag. Check it's formatted like a Python dictionary")
    return args

### Main prediction functions ###
def predict(traj, method, opts):
    """Predicts fraction of deuterium exchange for residues in the given
       trajectory, using the given method and dictionary of options.
 
       Usage: predict(traj, method, options)
       Returns: Object of desired method class, with completed HDX predictions"""

    # Switch for methods (add new ones here):
    methods = { 'radou' : Methods.Radou }

    result = methods[method.lower()](**opts)
    result.run(traj)
    return result

def full(trajlist, parm, stride, method, opts):
    """Loads all trajectories in the given list and performs HDX predictions.

       Usage: full(trajlist, parm, stride, method, options)
       Returns: Object of desired method class, with completed HDX predictions"""

    t = Functions.load_fulltraj(trajlist, parm=parm, stride=stride)
    return predict(t, method, opts)   

def chunks(trajlist, parm, stride, chunksize, method, opts):
    """Loads trajectories in the given list in chunks and performs HDX predictions.

       Usage: chunks(trajlist, parm, stride, chunksize, method, options)
       Returns: Object of desired method class, with completed HDX predictions"""

### Posibly rewrite this to sum result classes together into single object? '+=' = __iadd__(self,other)

    resultlist = []
    for t in trajlist:
        t_gen = Functions.load_trajchunks(t, parm=parm, stride=stride, chunk=chunksize)
        for t_chunk in t_gen:
            resultlist.append(predict(t_chunk, method, opts))

    return resultlist

def _update_options(opts, **updates):
    """Updates options dictionary with extra kwargs"""
    opts.update(updates)


### Main below here
if __name__ == '__main__':
    global args    
    args = parse()
    if args.method_options is not None:
        _update_options(args.method_options, logfile=args.logfile, \
                        segfile=args.segfile, outprefix=args.outprefix,\
                        times=args.times)
    else:
        args.method_options = {}
        _update_options(args.method_options, logfile=args.logfile, \
                        segfile=args.segfile, outprefix=args.outprefix,\
                        times=args.times)
    if args.chunks is not None:
        results = chunks(args.traj, args.parm, args.stride, args.chunks, args.method, args.method_options)
    else:
        results = full(args.traj, args.parm, args.stride, args.method, args.method_options)
        
