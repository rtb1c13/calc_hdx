#!/usr/bin/env python

# Author: Richard Bradshaw, richard.bradshaw@nih.gov

# Starting script for HDX analyses (could be updated into a library later?)

# Requirements - MDAnalysis, numpy

# Starts with script 1


import MDAnalysis as MD
import numpy as np
import argparse

# Flags for taking periodicity into account for contacts selections:
MD.core.flags['use_periodic_selections'] = True
MD.core.flags['use_KDTree_routines'] = False

# Functions
def setup_universe(parm, traj, pformat=None, tformat='DCD'):
    """Defines an MDAnalysis Universe object with the desired topology
       and coordinates.
       Usage: setup_universe(parm,traj,[pformat, tformat])
       pformat defaults to None which guesses parmfile formats
       based on their extension. tformat defaults to DCD for NAMD files. 
       If you're using a different format you'll have to set tformat yourself.
       Returns MDAnalysis Universe object"""
    return Universe(parm, traj, format=tformat, topology_format=pformat)


def list_prolines(univ, log="HDX_analysis.log"):
    """Creates a list of proline residues and appropriate resids

       Usage: list_prolines(univ, [log])
       Returns: Numpy array of [[Proline_ID, Proline_index]]"""
    try:
        prolist = [ r.resid for r in univ.PROT.residues if r.resname=='PRO' ]
        proidx = [ r.resindex for r in univ.PROT.residues if r.resname=='PRO' ]
    except AttributeError:
        with open(log, 'a') as f:
            f.write("No protein segment named 'PROT' found, \
                     searching all segments for prolines\n")
        prolist = [ r.resid for r in univ.residues if r.resname=='PRO' ]
        proidx = [ r.resindex for r in univ.residues if r.resname=='PRO' ]
    with open(log, 'a') as f:
        f.write("Prolines identified at resid:\n"+ \
                "%s\n" % ' '.join(str(i) for i in prolist))
    return np.asarray(zip(prolist, proidx))

def select_resids(univ, idxlist, protonly=True, invertgroup=None):
    """Returns atom selection of residue IDs in the supplied list, with options
       to restrict the selection to protein-only atoms (default) and/or 
       select all atoms NOT in a supplied Atom Group (inversion of selection, off
       by default)

       Usage: select_resids(univers, indexlist, [protonly, invertgroup])"""
    if inversion is not None:
        if protonly:
            return univ.select_atoms("protein and not group inversion", inversion=invertgroup)
        else:
            return univ.select_atoms("not group inversion", inversion=invertgroup)
    elif protonly:
        return univ.select_atoms("protein and resid %s" % ' '.join(str(r) for r in idxlist))
    else:
        return univ.select_atoms("resid %s" % ' '.join(str(r) for r in idxlist))

def calc_contacts(univ, selstr, cutoff):
    """Creates an atom selection of selected atoms within a specified cutoff"""
    univ.select_atoms(
