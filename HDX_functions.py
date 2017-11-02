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
                    "%s\n" % '\n'.join(atm2res(i) for i in traj.topology.select(atomselect))) 
        return traj.topology.select(atomselect)

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

def _calc_hbonds_bh(traj, HN, minfreq=0.0, cutoff=0.25, angle=120.0, **kwargs):
    """Calculates number of protein H-bonds for a particular atom index
       using the 'Baker-Hubbard' method. Default donor-acceptor distance < 0.25 nm
       + angle > 120 degrees.
       Reports all H-bonds (minimum freq=0.0) by default. For other kwargs
       see mdtraj.geometry.hbond._get_bond_triplets or ..._compute_bounded_geometry.
       
       Usage: _calc_hbonds_bh(traj, atom, [minfreq, cutoff, angle, **kwargs])
       Returns: n_frames length array of H-bond counts for desired atom"""

    # Atoms for H-bonds includes all O*, N* and single HN hydrogen
    
    c = traj.atom_slice(traj.topology.select("protein and (symbol O or symbol N) or index %s" % HN))

    # Call internal functions of md.baker_hubbard directly to return
    # distances & angles, otherwise only bond_triplets averaged across
    # a trajectory are returned
    bond_triplets = md.geometry.hbond._get_bond_triplets(c.topology, **kwargs)
    mask, distances, angles = md.geometry.hbond._compute_bounded_geometry(c, bond_triplets,\
                              0.25, [1, 2], [0, 1, 2], freq=minfreq, periodic=True)

    # can select distance/angle criteria here
    try:
        ang_rad = 2.0*np.pi / (360./angle)
    except ZeroDivisionError:
        angle = 360.0
        ang_rad = 2.0*np.pi / (360./angle)
        
    hbond_counts = np.sum(np.logical_and(distances < cutoff, angles > ang_rad), axis=1)
    return hbond_counts

def calc_hbonds(traj, method, donors, skip_first=True, **kwargs):
    """Calculates H-bond counts per frame for each atom in 'donors' array
       to each acceptor atom in the system. H-bonds can be defined using
       any one of the methods below.

       Available methods:
          'contacts' : Distance-based cutoff of 0.24 nm 
          'bh'       : Baker-Hubbard distance ( < 0.25 nm) and angle ( > 120 deg) cutoff

       Default cutoff/angle can be adjusted with the 'cutoff' and 'angle' kwargs

       Usage: calc_hbonds(traj, method=['contacts','bh'], donors, [**kwargs])
       Returns: n_donors * n_frames 2D array of H-bond counts per frame for all donors"""
    
    # Switch for H-bond methods
    methods = {
               'contacts' : _calc_hbonds_contacts,
               'bh' : _calc_hbonds_bh
              }

    if skip_first:
        donors = donors[1:] # Remove first atom/residue from list

    try:
        total_counts = np.zeros((len(donors),traj.n_frames))
    except TypeError:
        total_counts = np.zeros((1,traj.n_frames))
    for i, v in enumerate(donors):
        total_counts[i] = methods[method](traj, v, **kwargs)

    reslist = [ traj.topology.atom(a).residue.index for a in donors ]
#    hbonds = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), total_counts), axis=1) # Array of [[ Res idx, Contact count ]]

    return np.asarray(reslist), total_counts 
        
def calc_nh_contacts(traj, reslist, cutoff=0.65, skip_first=True, protonly=True):
    """Calculates contacts between each NH atom and the surrounding heavy atoms,
       excluding those in residues n-2 to n+2.

       By default contacts < 0.65 nm are calculated, and only protein-heavys,
       are included, but can include all heavys if desired. Also skips first 
       residue (N-terminus) in a residue list by default too.

       Usage: calc_nh_contacts(traj, reslist, [cutoff=0.65, skip_first=True, protein=True])
       Returns: n_res x n_frames 2D array of contacts per frame for each residue in reslist"""

    # Check if current atom is a heavy atom
    is_heavy = lambda _: traj.topology.atom(_).element.symbol is not 'H'

    if skip_first:
        reslist.pop(0) # Remove first atom/residue from list

    contact_count = np.zeros((len(reslist), traj.n_frames))
    for idx, res in enumerate(reslist):
        robj = traj.topology.residue(res)
        excl_idxs = range(robj.index - 2, robj.index + 3, 1) # Exclude n-2 to n+2 residues
        inv_atms = select_residxs(traj, excl_idxs, protonly=protonly, invert=True) # At this stage includes H + heavys
        heavys = inv_atms[ np.array( [ is_heavy(i) for i in inv_atms ] ) ] # Filter out non-heavys
        
        contact_count[idx] = calc_contacts(traj, robj.atom('N').index, heavys, cutoff=cutoff)

#    contacts = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), contact_count), axis=1) # Array of [[ Res idx, Contact count ]]
    return np.asarray(reslist), contact_count 

def PF(traj, hbond_method='contacts', save_contacts=False, **kwargs):
    """Calculates Radou et al. protection factors for a provided trajectory.
       Empirical scaling factors of Nh * 2.0 and Nc * 0.35, as per Radou et al.
       H-bonds can be calculated using either the 'contacts' definition (all
       acceptors within 0.24 nm by default) or the Baker-Hubbard distance +
       angle definition (0.25 nm / 120 deg by default). Printout of temporary
       files containing by-residue contacts can be enabled with save_contacts.

       All proline residues and the N-terminal residue are skipped. See 
       calc_hbonds and calc_nh_contacts for optional kwargs.       

       Usage: PF(traj, [ hbond_method=['contacts','bh'], save_contacts=False, **kwargs ])
       
       Returns: (array of residue indices, array of mean protection factors)"""
    # Setup residue/atom lists        
    hn_atms = extract_HN(traj)
    prolines = list_prolines(traj)
    # All protein residues except prolines
    reslist = [ r.index for r in traj.topology.residues if r.is_protein and r.index not in prolines[:,1] ]

    # Calc Nc/Nh
    hres, hbonds = calc_hbonds(traj, hbond_method, hn_atms, **kwargs)
    cres, contacts = calc_nh_contacts(traj, reslist, **kwargs)

    if not np.array_equal(hres, cres):
        raise HDX_Error("The residues analysed for Nc and Nh appear to be different. Check your inputs!")

    # Option to save outputs
    if save_contacts:
        for i, residx in enumerate(hres):
            np.savetxt("Hbonds_%d.tmp" % (residx+1), hbonds[i], fmt='%d') # Use residue indices internally, print out IDs
        for i, residx in enumerate(cres):
            np.savetxt("Contacts_%d.tmp" % (residx+1), contacts[i], fmt='%d') # Use residue indices internally, print out IDs
    # Calc PF with phenomenological equation
    hbonds *= 2        # Beta parameter 1
    contacts *= 0.35   # Beta parameter 2

    pf = np.exp(hbonds + contacts)
    pf = np.mean(pf, axis=1)

    # Save PFs to separate log file
    np.savetxt("Protection_factors.tmp", np.stack((hres+1, pf), axis=1), \
               fmt=['%7d','%18.8f'], header="ResID  Protection factor") # Use residue indices internally, print out IDs

    return hres, pf

# Kints? First check cis/trans pro

#atoms, angles = md.compute_omega(traj)

def pro_omega_indices(traj, prolines):
    """Calculates omega dihedrals (CA-C-N-CA) for all proline
       residues in a given prolines array from list_prolines.

       Usage: pro_omega_indices(traj, prolines)
       Returns: (atom_indices, w_angles_by_frame)"""

    atom_names = ['CA', 'C', 'N', 'CA']
    offset = np.asarray([-1, -1, 0, 0])

    w_indices = np.zeros((len(prolines),4))
    for i, residx in enumerate(prolines[:,1]):
        curr_offset = offset + residx
        curr_atm_indices = []
        for x in zip(curr_offset, atom_names):
            # Cycle through previous CA/C, current N, CA
            curr_atm_indices.append(traj.topology.residue(x[0]).atom(x[1]).index)
        w_indices[i] = curr_atm_indices

    return w_indices, md.compute_dihedrals(traj, w_indices)
        
# Assignments for intrinsic rate adjustments:
# 1) cis/trans prolines
# 2) disulfides
# 3) His protonation states
# 4) N/C termini
       
def assign_cis_proline(traj, log="HDX_analysis.log"):
    """Assigns cis-proline residues"""

    prolines = list_prolines(traj)
    outidxs, outangs = pro_omega_indices(traj, prolines)
    for i, proidx in enumerate(prolines[:,1]):
        traj.topology.residue(proidx).cis_byframe = np.logical_and(outangs < np.pi/2, outangs > -1*np.pi/2)[:,i]
        if np.max(traj.topology.residue(proidx).cis_byframe) > 0:
            with open(log, 'a') as f:
               f.write("Cis-proline found at frame %d for residue %s!\n" % (np.argmax(traj.topology.residue(proidx).cis_byframe) + 1, traj.topology.residue(proidx).resSeq))


def assign_disulfide(traj, log="HDX_analysis.log"):
    """Assigns residues involved in disulfide bridges"""

    # This selection syntax & assignment is untested
    sg = traj.topology.select('name SG and resname CYS or resname CYX')
    try:
        sg_coords = traj.atom_slice(sg).xyz
    except IndexError: # Catch empty sg
        with open(log, 'a') as f:
            f.write("No cysteines found in topology.\n")
        return
    traj.topology.create_standard_bonds()
    traj.topology.create_disulfide_bonds(sg_coords)
    # Assign disulfides (identical for each frame)
    for b in traj.topology._bonds:
        if all(i.element.symbol == 'S' for i in b):
            b[0].residue.disulf = np.ones(traj.n_frames, dtype=bool)
            b[1].residue.disulf = np.ones(traj.n_frames, dtype=bool)
            with open(log, 'a') as f:
                f.write("Disulfide found for residues %s - %s\n" \
                         % (b[0].residue, b[1].residue))
    
def assign_his_protonation(traj, log="HDX_analysis.log"):
    """Assigns protonation state to HIS residues"""

    hisidx = [ r.index for r in traj.topology.residues if r.code == 'H' ]
    for i in hisidx:
        names = [ a.name for a in traj.topology.residue(i).atoms ]
        if all(n in names for n in ['HD1','HE2']): # Atom names for doubly protonated His (Hip)
            traj.topology.residue(i).HIP = np.ones(traj.n_frames, dtype=bool)
            with open(log, 'a') as f:
                f.write("Protonated His assigned for residue %d\n" % traj.topology.residue(i).resSeq)
#        else:
#            traj.topology.residue(i).HIP = np.zeros(traj.n_frames, dtype=bool)

def assign_termini(traj, log="HDX_analysis.log"):
    """Assigns flags to N and C terminal residues"""

    for c in traj.topology.chains:
        c.residue(0).nterm = np.ones(traj.n_frames, dtype=bool)
        c.residue(-1).cterm = np.ones(traj.n_frames, dtype=bool)
        with open(log, 'a') as f:
            f.write("N-terminus identified at: %s\nC-terminus identified at: %s\n" \
                     % (c.residue(0), c.residue(-1)))
                     
# Dict for intrinsic rate calculations

# Each key is a residue name, each value is [ lgAL, lgAR, lgBL, lgBR ]
### THIS IS A DIFFERENT ORDER TO TABLE 2 IN THE REFERENCE BELOW ###
# Values from Bai et al., Proteins, 1993, 17, 75-86


# Note that Englander has adjustments to Glu/Asp rates in spreadsheet
# on his website, based on Mori et al., Proteins, 1997, 28, 325-332

rate_adjs = { 'ALA': [ 0.00, 0.00, 0.00, 0.00 ],
              'ARG': [ -0.59, -0.32, 0.08, 0.22 ],         
              'ASN': [ -0.58, -0.13, 0.49, 0.32 ],         
              'ASP': [ 0.90, 0.58, -0.30, -0.18 ],         
              'ASH': [ -0.90, -0.12, 0.69, 0.60 ], # Protonated ASP        
              'CYS': [ -0.54, -0.46, 0.62, 0.55 ],         
              'CYS2': [ -0.74, -0.58, 0.55, 0.46 ], # Disulfide         
              'GLY': [ -0.22, 0.22, 0.27, 0.17 ],         
              'GLN': [ -0.47, -0.27, 0.06, 0.20 ],         
              'GLU': [ -0.90, 0.31, -0.51, -0.15 ],         
              'GLH': [ -0.60, -0.27, 0.24, 0.39 ], # Protonated GLU        
              'HIS': [ 0.00, 0.00, -0.10, 0.14 ],  # Acid rates are N/D, 
                                                   # but at pH where His is deprotonated,
                                                   # errors will be negligible       
              'HIP': [ -0.80, -0.51, 0.80, 0.83 ],         
              'ILE': [ -0.91, -0.59, -0.73, -0.23 ],         
              'LEU': [ -0.57, -0.13, -0.58, -0.21 ],         
              'LYS': [ -0.56, -0.29, -0.04, 0.12 ],         
              'MET': [ -0.64, -0.28, -0.01, 0.11 ],         
              'PHE': [ -0.52, -0.43, -0.24, 0.06 ],         
              'PRO': [ 0.00, -0.19, 0.00, -0.24 ], # Trans PRO        
              'PROC': [ 0.00, -0.85, 0.00, 0.60 ], # Cis PRO        
              'SER': [ -0.44, -0.39, 0.37, 0.30 ],         
              'THR': [ -0.79, -0.47, -0.07, 0.20 ],         
              'TRP': [ -0.40, -0.44, -0.41, -0.11 ],         
              'TYR': [ -0.41, -0.37, -0.27, 0.05 ],         
              'VAL': [ -0.74, -0.30, -0.70, -0.14 ],         
              'NT': [ 0.00, -1.32, 0.00, 1.62 ], # N-term NH3+        
              'CT': [ 0.96, 0.00, -1.80, 0.00 ], # C-term COO-        
              'CTH': [ 0.05, 0.00, 0.00, 0.00 ], # C-term COOH
            }
# Adjust ordering so value is [ lgAL, lgBL, lgAR, lgBR ] for kint analysis 
### THIS IS A DIFFERENT ORDER TO TABLE 2 IN THE REFERENCE ABOVE ###
_rate_adjs = rate_adjs.copy()
for i in _rate_adjs.values(): 
    i[1], i[2] = i[2], i[1]     # Swap elements


# Helper function to turn sequence-specific rate adjustments to intrinsic acid/base/water rates
def _adj_to_rates(rate_adjs, lgkAref=2.04, lgkBref=10.36, lgkWref=-1.5, \
                  EaA=14., EaB=17., EaW=19., R=0.001987, Tref=293, \
                  Texp=298, pKD=14.87, pD=7.4):
    """Calculates intrinsic rates for a given set of rate adjustments
       [ log(AL), log(BL), log(AR), log(BR) ] taken from Bai et al."""

    # Calc reference rates at experimental temperature
    # / np.log(10) = conversion from ln to log10
    lgkAexp = lgkAref - (EaA / np.log(10) / R) * (1./Texp - 1./Tref)
    lgkBexp = lgkBref - (EaB / np.log(10) / R) * (1./Texp - 1./Tref)
    lgkWexp = lgkWref - (EaW / np.log(10) / R) * (1./Texp - 1./Tref)

    # Calc log(kA||kB||kW)
    lgkA = lgkAexp + rate_adjs[0] + rate_adjs[2] - pD 
    lgkB = lgkBexp + rate_adjs[1] + rate_adjs[3] - pKD + pD  
    lgkW = lgkWexp + rate_adjs[1] + rate_adjs[3] 

    kint = 10**lgkA + 10**lgkB + 10**lgkW
    #print(lgkAexp, lgkBexp, lgkWexp, 10**lgkA, 10**lgkB, 10**lgkW)
    return kint

# Intrinsic rate calc:
def kint(traj, reslist, log="HDX_analysis.log", **kwargs):
    """Function for calculating intrinsic rates of residues
       in a given topology
       
       Intrinsic exchange rates k_int are computed using equations below.
       k_int = k_A + k_B + k_W
       lgk_A = lgk_A,ref + lgA_L + lgA_R - pD
       lgk_B = lgk_B,ref + lgB_L + lgB_R - pOD
             = lgk_B,ref + lgB_L + lgB_R - pK_D + pOD
       lgk_W = lgk_W,ref + lgB_L + lgB_R"""
    # Create tmp copy of topology to adjust residue names
    tmptop = traj.topology.copy()
    
    kints = np.zeros(len(reslist))

    # Adjust residue names for: Cis-Pro, HIP, cystine bridges, GLH/ASH
    reslist = np.insert(reslist,0,reslist[0]-1) # Insert 'prev' residue for first index
    for i in reslist:
        curr = tmptop.residue(i)
        try:
            if np.max(curr.cis_byframe): # If cis-proline is true for any frame
                curr.name = 'PROC'
                continue
        except AttributeError:
            pass
        try:
            if np.max(curr.HIP): # If His+ is true for any frame
                curr.name = 'HIP'
                continue
        except AttributeError:
            pass
        try:
            if np.max(curr.disulf): # If Cys-Cys is true for any frame
                curr.name = 'CYS2'
                continue
        except AttributeError:
            pass
        if curr.name == 'GLU': # If Glu has a prodonated carboxylate
            try:
                curr.atom('HE2')
                curr.name = 'GLH'
                continue
            except KeyError:
                pass
        if curr.name == 'ASP': # If Asp has a protonated carboxylate
            try:
                curr.atom('HD2')
                curr.name = 'ASH'
                continue
            except KeyError:
                pass

    # Assign N/C termini
    for chain in tmptop.chains:
        chain.residue(0).name = 'NT'
        # Doesn't appead to be a standard name for COOH hydrogen - guess based on no. of bonds!
        if chain.residue(-1).atom('O').n_bonds > 1 or chain.residue(-1).atom('OXT').n_bonds > 1:
            chain.residue(-1).name = 'CTH'
            print("It looks like you have a neutral C-terminus (COOH) at residue %s?" % chain.residue(-1)) 
        else:
            chain.residue(-1).name = 'CT'
        
    reslist = np.delete(reslist, 0) # Remove i-1 residue we inserted above
    for i, r in enumerate(reslist):
        curr = tmptop.residue(r)
        prev = tmptop.residue(r-1)
        # check for cispro
        if prev.name == 'PROC':
            with open(log, 'a') as f:
                f.write("Performing rate calculation on a by-frame basis for residue %s" % prev)
            adj_cis, adj_trans = _rate_adjs[curr.name][0:2], _rate_adjs[curr.name][0:2] 
            adj_cis.extend(_rate_adjs['PROC'][2:4])
            adj_trans.extend(_rate_adjs['PRO'][2:4])
            kint_cis = _adj_to_rates(adj_cis, **kwargs)
            kint_trans = _adj_to_rates(adj_trans, **kwargs)
            kint_byframe = np.where(prev.cis_byframe, kint_cis, kint_trans)
            kints[i] = np.mean(kint_byframe) # Overall intrinsic rate is adjusted by mean population of cis-pro
        else:
            curr_adjs = _rate_adjs[curr.name][0:2]
            curr_adjs.extend(_rate_adjs[prev.name][2:4])
            kints[i] = _adj_to_rates(curr_adjs, **kwargs)

    # Save PFs to separate log file
    np.savetxt("Intrinsic_rates.tmp", np.stack((reslist+1, kints), axis=1), \
               fmt=['%7d','%18.8f'], header="ResID  Intrinsic rate ") # Use residue indices internally, print out IDs

    return kints, reslist

# main() below here
#if __name__ == '__main__-:
     
