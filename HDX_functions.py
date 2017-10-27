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
       atoms to keep) and stride (integer of every nth frame to keep).

       'standard_names=False' may also be useful for PDB topologies,
       otherwise amide H might be renamed from the atom names provided
       to the standard PDB identifiers (e.g. 'H', 'H2', 'H3' for the 
       terminal NH3 group).
       
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
       
       'standard_names=False' may also be useful for PDB topologies,
       otherwise amide H might be renamed from the atom names provided
       to the standard PDB identifiers (e.g. 'H', 'H2', 'H3' for the 
       terminal NH3 group).
   
       Returns a generator object with trajectory iterations."""
    return md.iterload(traj, top=parm, **kwargs)



def list_prolines(traj, log="HDX_analysis.log"):
    """Creates a list of proline residues and appropriate resids.
       Resids are output to HDX_analysis.log file by default.

       Usage: list_prolines(traj, [log])
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
       Periodicity is included in MDtraj function by default.
       Usage: calc_contacts(traj, qidx, cidx, [cutoff=0.65]).

       Qidx and cidx are the atom index lists to search for contacts from
       and to respectively (e.g. from amide NH to all heavy atoms).

       Returns count of contacts for each frame in supplied trajectory."""

    try:
        byframe_ctacts = md.compute_neighbors(traj, cutoff, qidx, haystack_indices=cidx)
    except TypeError:
        print("Now calculating contacts to single atom, idx %d" % qidx)
        qidx = np.array([qidx])
        byframe_ctacts = md.compute_neighbors(traj, cutoff, qidx, haystack_indices=cidx)
    return map(lambda x: len(x), byframe_ctacts)


# These should go in a class for residues that calcs pfactor, kint, neighbouring residues and the rest
def extract_HN(traj, prolines=None, atomselect="(name H or name HN)", log="HDX_analysis.log"):
    """Returns a list of backbone amide H atom indices, suitable
       for use with 'calc_contacts'. Optionally takes an array of 
       resids/indices to skip (normally prolines) and by default returns
       atom indices matching 'name H and backbone'
       
       Usage: extract_NH(traj, [prolines, atomselect])"""

    # Combine res name & ID to concatenated identifier
    atm2res = lambda _: traj.topology.atom(_).residue.name + str(traj.topology.atom(_).residue.resSeq)

    if prolines is not None:
        # Syntax = "... and not (residue 1 or residue 2 or residue 3 ... )"
        atomselect += " and not (residue %s" % ' or residue '.join(str(_) for _ in prolines[:,0]) + ")"
        with open(log, 'a') as f:
            f.write("Extracted HN from resids:\n"+ \
                    "%s\n" % '\n'.join(atm2res(i) for i in traj.topology.select(atomselect)))
        return traj.topology.select(atomselect)
    else:
        with open(log, 'a') as f:
            f.write("Extracted HN from resids:\n"+ \
                    "%s\n" % ' '.join(atm2res(i) for i in traj.topology.select(atomselect)))
        return traj.topology.select(atomselect)

   # Next: take each HN and calc hbonds?
   # calc_contacts can be used for normal distance-based selections (2.4A H - heavy)
def _calc_hbonds_contacts(traj, HN, cutoff=0.24, **kwargs):
    """Calculates number of protein H-bonds for a particular atom index
       using the 'contacts' method. Cutoff = 0.24 nm by default. Bonds to
       all O* or N* evaluated.
       
       Usage: _calc_hbonds_contacts(t, atom, [cutoff])"""

    # Get N index in same residue as current HN atom
    getN4H = lambda _: traj.topology.atom(_).residue.atom('N').index
    c = traj.topology.select("protein and (symbol O or symbol N) and not index %s" % getN4H(HN))

    hbond_counts = calc_contacts(traj, HN, c, cutoff=cutoff)
    return hbond_counts

def _calc_hbonds_bh(traj, HN, minfreq=0.0, **kwargs):
    """Calculates number of protein H-bonds for a particular atom index
       using the 'Baker-Hubbard' method (donor-acceptor distance < 0.25 nm
       + angle > 120 degrees).
       Reports all H-bonds (minimum freq=0.0) by default. Other kwargs:
       exclude_water=True and periodic=True.
       
       Usage: _calc_hbonds_bh(t, atom, [minfreq, **kwargs])"""

    # Atoms for H-bonds includes all O*, N* and single HN hydrogen
    
    c = traj.atom_slice(traj.topology.select("protein and (symbol O or symbol N) or index %s" % HN))
    bond_triplets = md.geometry.hbond._get_bond_triplets(c.topology,\
                    exclude_water=True, sidechain_only=False)
    def _compute_bounded_geometry_RTB(traj, triplets, distance_cutoff, distance_indices,
                                      angle_indices, freq=0.0, periodic=True):
        """
        Returns a tuple include (1) the mask for triplets that fulfill the distance
        criteria frequently enough, (2) the actual distances calculated, and (3) the
        angles between the triplets specified by angle_indices.
        """
        # First we calculate the requested distances
        distances = md.compute_distances(traj, triplets[:, distance_indices], periodic=periodic)

        # Calculate angles using the law of cosines
        abc_pairs = zip(angle_indices, angle_indices[1:] + angle_indices[:1])
        abc_distances = []

        # Calculate distances (if necessary)
        for abc_pair in abc_pairs:
            if set(abc_pair) == set(distance_indices):
                abc_distances.append(distances)
            else:
                abc_distances.append(md.compute_distances(traj, triplets[:, abc_pair],
                    periodic=periodic))

        # Law of cosines calculation
        a, b, c = abc_distances
        cosines = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
        np.clip(cosines, -1, 1, out=cosines) # avoid NaN error
        angles = np.arccos(cosines)

        return distances, angles
    distances, angles = _compute_bounded_geometry_RTB(c, bond_triplets,\
                              0.25, [1, 2], [0, 1, 2], freq=0.0, periodic=True)
    hbond_counts = np.sum(np.logical_and(distances < 0.25, angles > (2.0 * np.pi / 3)),axis=1)
    return hbond_counts

def calc_hbonds(traj, method, donors, **kwargs):
    """Calculates H-bond counts per frame for each atom in 'donors' array
       to each acceptor atom in the system. H-bonds can be defined using
       any one of the methods below.

       Available methods:
          'contacts' : Distance-based cutoff of 0.24 nm
          'bh'       : Baker-Hubbard distance ( < 0.25 nm) and angle ( > 120 deg) cutoff
       
       Usage: calc_hbonds(traj, method=['contacts','bh'], donors)"""
    
    # Switch for H-bond methods
    methods = {
               'contacts' : _calc_hbonds_contacts,
               'bh' : _calc_hbonds_bh
              }
    try:
        total_counts = np.zeros((len(donors),traj.n_frames))
    except TypeError:
        total_counts = np.zeros((1,traj.n_frames))
    for i, v in enumerate(donors):
        total_counts[i] = methods[method](traj, v, **kwargs)
    return total_counts
        
          


   # Then: Take each NH and calc contacts?

   # Then: PF calculation function?







 
