#!/usr/bin/env python

# Transferable functions for HDX analysis

import mdtraj as md
import numpy as np


# Exception for HDX
class HDX_Error(Exception):
    """Exception in HDX module"""

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


def select_residxs(traj, reslist, protonly=True, invert=False):
    """Returns atom indices of atoms belonging to residues (0-indexed)
       in the supplied list,

       Options to restrict the selection to protein-only atoms
       (default) and/or select all atoms NOT in residues in the supplied list.
       (inversion of selection, off by default)

       Usage: select_resids(traj, reslist, [protonly, invert])
       Returns: Numpy array of selected atom indices"""

    # The topology.select syntax is more restrictive than MDAnalysis here
    # - use list comprehensions instead
    if invert:
        if protonly:
            return np.asarray([ atom.index for atom in traj.topology.atoms if (atom.residue.is_protein and atom.residue.index not in reslist) ])
        else:
            return np.asarray([ atom.index for atom in traj.topology.atoms if (atom.residue.index not in reslist) ])
    elif protonly:
        return np.asarray([ atom.index for atom in traj.topology.atoms if (atom.residue.is_protein and atom.residue.index in reslist) ])
    else:
        return np.asarray([ atom.index for atom in traj.topology.atoms if (atom.residue.index in reslist) ])


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
                    "%s\n" % '\n'.join(atm2res(i) for i in traj.topology.select(atomselect))) 
        return traj.topology.select(atomselect)


