#!/usr/bin/env python

# Author: Richard Bradshaw, richard.bradshaw@nih.gov

# Starting script for HDX analyses (could be updated into a library later?)

# Requirements - MDTraj, numpy, argparse

# Starts with script 1

# For compatibility with Patrick's python 3
from __future__ import print_function
from __future__ import division
#
import mdtraj as md
import numpy as np
import argparse



# Functions
def load_fulltraj(traj, parm, **kwargs):
    """Loads an MDtraj trajectory object with the desired topology
       and coordinates.
       Usage: setup_universe(parm,traj,[**kwargs])
       Standard kwargs include atom_indices (an array of 0-indexed
       atoms to keep) and stride (integer of every nth frame to keep) 
       Returns a complete trajectory, which may be memory intensive.

       See also load_trajchunks for an iterative load of large trajectories """
    return md.load(traj, top=parm, **kwargs)

 def load_trajchunks(traj, parm, **kwargs):
    """Loads a file into a generator of MDtraj trajectory chunks.
       Useful for large/memory intensive trajectory files
       Usage: load_trajchunks(traj, parm, [**kwargs])
       Standard kwargs include chunk (size of the trajectory chunks
       to load per iteration), skip (n frames to skip at the start of
       the trajectory), atom_indices (an array of 0-indexed
       atoms to keep) and stride (integer of every nth frame to keep) 
   
       Returns a generator object with trajectory iterations."""
    return md.iterload(traj, top=parm, **kwargs)



def list_prolines(univ, log="HDX_analysis.log"):
    """Creates a list of proline residues and appropriate resids

       Usage: list_prolines(univ, [log])
       Returns: Numpy array of [[Proline_ID, Proline_index]]"""
    prolist = [ r.resSeq for r in traj.topology.residues if r.name=='PRO' ]
    proidx = [ r.index for r in traj.topology.residues if r.name=='PRO' ]
    with open(log, 'a') as f:
        f.write("Prolines identified at resid:\n"+ \
                "%s\n" % ' '.join(str(i) for i in prolist))
    return np.asarray(zip(prolist, proidx))

def select_resids(traj, idxlist, protonly=True, invert=False):
    """Returns atom indices of residues (0-indexed) in the supplied list,
       with options to restrict the selection to protein-only atoms
       (default) and/or select all atoms NOT in a supplied list.
       (inversion of selection, off by default)

       Usage: select_resids(traj, residxlist, [protonly, invert])"""

    # The topology.select syntax is more restrictive than MDAnalysis here
    # - use list comprehensions instead
    if invert:
        if protonly:
            return np.asarray([ atom.index for atom in t.topology.atoms if (atom.residue.is_protein and atom.residue.index not in l) ])
        else:
            return np.asarray([ atom.index for atom in t.topology.atoms if (atom.residue.index not in l) ])
    elif protonly:
        return np.asarray([ atom.index for atom in t.topology.atoms if (atom.residue.is_protein and atom.residue.index in l) ])
    else:
        return np.asarray([ atom.index for atom in t.topology.atoms if (atom.residue.index in l) ])

def calc_contacts(traj, qidx, cidx, cutoff=0.65):
    """Calculates contacts between 'query' and 'contact' atom selections
       within a specified cutoff (default = 0.65, for coordinates in nm).
       Usage: calc_contacts(traj, qidx, cidx, [cutoff=0.65]).

       Qidx and cidx are the atom index lists to search for contacts from
       and to respectively (e.g. from amide NH to all heavy atoms)."""
    
