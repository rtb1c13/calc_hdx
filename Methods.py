#!/usr/bin/env python

# Class for HDX trajectories, inherited from MDTraj
from __future__ import print_function
from __future__ import division
# 
import mdtraj as md
import numpy as np
import os, glob, itertools
import Functions


class Radou():
    """Class for Radou-style analysis. Initialises with a dictionary of default
       parameters for analysis, accessible as Radou.params

       Default parameters can either be updated directly in the Radou.params
       dictionary or by supplying a extra parameters as kwargs during
       initialisation, e.g.: Radou(cut_nc=1.0) or Radou(**param_dict)

       Run a by-residue deuterated fraction prediction with these parameters
       using the Radou.run method."""

    def __init__(self, **extra_params):
        """Initialises parameters for Radou-style analysis.
           See self.params for default values"""
        # Initialise main parameters with defaults
        self.params = { 'hbond_method' : 'contacts',
                        'protonly' : True,
                        'cut_Nc' : 0.65,
                        'cut_Nh' : 0.24,
                        'bh_dist' : 0.25,
                        'bh_ang' : 120.0,
                        'save_contacts' : False,
                        'skip_first' : True,
                        'betac' : 0.35,
                        'betah' : 2.0,
                        'kint_adjs' : None,
                        'kint_params' : None,
                        'times' : [ 0.167, 1.0, 10.0, 120.0],
                        'segfile' : "cropped_seg.list",
                        'logfile' : "HDX_analysis.log",
                        'outprefix' : '' }
        self.params.update(extra_params) # Update main parameter set from kwargs

        # Default rate adjustments for adjacent amino acids
        # Each key is a residue name, each value is [ lgAL, lgAR, lgBL, lgBR ]
        # Values from Bai et al., Proteins, 1993, 17, 75-86

        # Note that Englander has adjustments to Glu/Asp rates in spreadsheet
        # on his website, based on Mori et al., Proteins, 1997, 28, 325-332
        # These are NOT included here by default
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
                      'CTH': [ 0.05, 0.00, 0.00, 0.00 ], } # C-term COOH

        # Add updates from kint_adjs dictionary if it's given as a kwarg
        if self.params['kint_adjs'] is not None:
            rate_adjs.update(self.params['kint_adjs'])
        self.params['kint_adjs'] = rate_adjs

        # Adjust ordering so value is [ lgAL, lgBL, lgAR, lgBR ] for kint analysis 
        ### THIS IS A DIFFERENT ORDER TO TABLE 2 IN THE REFERENCE ABOVE ###
        _reordered_rate_adjs = { k : v[:] for k, v in self.params['kint_adjs'].items() } # Deep copy
        for i in _reordered_rate_adjs.values():
            i[1], i[2] = i[2], i[1]
        self.params['_reordered_kint_adjs'] = _reordered_rate_adjs
 
        # Default parameters for ka/kb/kw estimations
        # Values from Bai et al., Proteins, 1993, 17, 75-86
        rate_params = { 'lgkAref' : 2.04,
                        'lgkBref' : 10.36,
                        'lgkWref' : -1.5,
                        'EaA' : 14.,
                        'EaB' : 17.,
                        'EaW' : 19.,
                        'R' : 0.001987,
                        'Tref' : 293,
                        'Texp' : 298,
                        'pKD' : 14.87,
                        'pD' : 7.4 }

        # Add updates from kint_adjs dictionary if it's given as a kwarg
        if self.params['kint_params'] is not None:
            rate_params.update(self.params['kint_params'])
        self.params['kint_params'] = rate_params

    def __str__(self):
        """Print the method name"""
        return 'Radou'

    def __iadd__(self, other):
        """Sum results in other method object to this one, weighted by number of frames in each"""
        if isinstance(other, Radou):
            try:
                if np.array_equal(self.rates, other.rates):
                    self.pfs = (self.n_frames * self.pfs) + (other.n_frames * other.pfs)
                    self.n_frames += other.n_frames
                    self.pfs /= self.n_frames
                    self.resfracs = self.dfrac(write=False)
                    return self
                else:
                    raise Functions.HDX_Error("Cannot sum two Radou objects with different intrinsic rates.")
            except AttributeError:
                return self
        else:
            return self

    def calc_contacts(self, qidx, cidx, cutoff):
        """Calculates contacts between 'query' and 'contact' atom selections
           within a specified cutoff (in nm).
           Periodicity is included in MDtraj function by default.
           Usage: calc_contacts(qidx, cidx, cutoff).

           Qidx and cidx are the atom index lists to search for contacts from
           and to respectively (e.g. from amide NH to all heavy atoms).

           Returns count of contacts for each frame in trajectory Radou.t."""

        try:
            byframe_ctacts = md.compute_neighbors(self.t, cutoff, qidx, haystack_indices=cidx)
        except TypeError:
#            print("Now calculating contacts to single atom, idx %d" % qidx)
            qidx = np.array([qidx])
            byframe_ctacts = md.compute_neighbors(self.t, cutoff, qidx, haystack_indices=cidx)
        return map(lambda x: len(x), byframe_ctacts)

    def _calc_hbonds_contacts(self, HN):
        """Calculates number of protein H-bonds for a particular atom index
           using the 'contacts' method. Bonds to all protein O* or N* evaluated
           by default, optionally all non-protein too (including waters) if 
           Radou.params['protonly'] is False.
       
           Usage: _calc_hbonds_contacts(atom)"""

        # Get N index in same residue as current HN atom
        getN4H = lambda _: self.top.atom(_).residue.atom('N').index
        if self.params['protonly']:
            c = self.top.select("protein and (symbol O or symbol N) and not index %s" % getN4H(HN))
        else:
            c = self.top.select("(symbol O or symbol N) and not index %s" % getN4H(HN))

        hbond_counts = self.calc_contacts(HN, c, self.params['cut_Nh'])
        return hbond_counts

    def _calc_hbonds_bh(self, HN, minfreq=0.0):
        """Calculates number of protein H-bonds for a particular atom index
           using the 'Baker-Hubbard' method. Default donor-acceptor distance < 0.25 nm
           + angle > 120 degrees in Radou.params.
           Reports all H-bonds (minimum freq=0.0) by default. Bonds to all protein 
           O* or N* evaluated by default, optionally all non-protein too 
           (including waters) if Radou.params['protonly'] is False.
       
           Usage: _calc_hbonds_bh(atom, [minfreq])
           Returns: n_frames length array of H-bond counts for desired atom"""

        # Atoms for H-bonds includes protein or all O*, N* and single HN hydrogen

        if self.params['protonly']:
            c = self.t.atom_slice(self.top.select("protein and (symbol O or symbol N) or index %s" % HN))
        else:
            c = self.t.atom_slice(self.top.select("all (symbol O or symbol N) or index %s" % HN))

        # Call internal functions of md.baker_hubbard directly to return
        # distances & angles, otherwise only bond_triplets averaged across
        # a trajectory are returned
        bond_triplets = md.geometry.hbond._get_bond_triplets(c.topology, exclude_water=self.params['protonly'])
        mask, distances, angles = md.geometry.hbond._compute_bounded_geometry(c, bond_triplets,\
                                  self.params['bh_dist'], [1, 2], [0, 1, 2], freq=minfreq, periodic=True)

        # can select distance/angle criteria here
        try:
            ang_rad = 2.0*np.pi / (360./self.params['bh_ang'])
        except ZeroDivisionError:
            
            self.params['bh_ang'] = 360.0
            ang_rad = 2.0*np.pi / (360./self.params['bh_ang'])

        hbond_counts = np.sum(np.logical_and(distances < self.params['bh_dist'], angles > ang_rad), axis=1)
        return hbond_counts

    def calc_hbonds(self, donors):
        """Calculates H-bond counts per frame for each atom in 'donors' array
           to each acceptor atom in the system. H-bonds can be defined using
           any one of the methods below, selected with Radou.params['hbond_method']
    
           Available methods:
              'contacts' : Distance-based cutoff of 0.24 nm 
              'bh'       : Baker-Hubbard distance ( < 0.25 nm) and angle ( > 120 deg) cutoff

           Default cutoff/angle can be adjusted with entries 'cut_Nh'/'bh_dist'/
           'bh_ang'in Radou.params.

           Usage: calc_hbonds(donors)
           Returns: n_donors * n_frames 2D array of H-bond counts per frame for all donors"""

    # Switch for H-bond methods
        methods = {
                   'contacts' : self._calc_hbonds_contacts,
                   'bh' : self._calc_hbonds_bh
                  }

        if self.params['skip_first']:
            donors = donors[1:] # Remove first atom/residue from list

        try:
            total_counts = np.zeros((len(donors), self.t.n_frames))
        except TypeError:
            total_counts = np.zeros((1, self.t.n_frames))
        for i, v in enumerate(donors):
            total_counts[i] = methods[self.params['hbond_method']](v)

        reslist = [ self.top.atom(a).residue.index for a in donors ]
#        hbonds = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), total_counts), axis=1) # Array of [[ Res idx, Contact count ]]

        return np.asarray(reslist), total_counts

    def calc_nh_contacts(self, reslist):
        """Calculates contacts between each NH atom and the surrounding heavy atoms,
           excluding those in residues n-2 to n+2.
    
           By Radou.params default contacts < 0.65 nm are calculated, and only
           protein-heavys, are included, but can include all heavys if desired.
           Also skips first residue (N-terminus) in a residue list by default too
           - see Radou.params['protonly'] and Radou.params['skip_first']

           Usage: calc_nh_contacts(reslist)
           Returns: (reslist, n_res x n_frames 2D array of contacts per frame for each residue)"""

        # Check if current atom is a heavy atom
        is_heavy = lambda _: self.top.atom(_).element.symbol is not 'H'

        if self.params['skip_first']:
            reslist.pop(0) # Remove first atom/residue from list

        contact_count = np.zeros((len(reslist), self.t.n_frames))
        for idx, res in enumerate(reslist):
            robj = self.top.residue(res)
            excl_idxs = range(robj.index - 2, robj.index + 3, 1) # Exclude n-2 to n+2 residues

            inv_atms = Functions.select_residxs(self.t, excl_idxs, protonly=self.params['protonly'], invert=True) # At this stage includes H + heavys
            heavys = inv_atms[ np.array( [ is_heavy(i) for i in inv_atms ] ) ] # Filter out non-heavys

            contact_count[idx] = self.calc_contacts(robj.atom('N').index, heavys, cutoff=self.params['cut_Nc'])

#        contacts = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), contact_count), axis=1) # Array of [[ Res idx, Contact count ]]
        return np.asarray(reslist), contact_count

    def PF(self):
        """Calculates Radou et al. protection factors for a provided trajectory.
           Empirical scaling factors of Nh * betah and Nc * betac taken from 
           Radou.params (2.0 & 0.35 respectively by default).
           H-bonds can be calculated using either the 'contacts' definition or
           the Baker-Hubbard distance + angle definition. Printout of temporary
           files containing by-residue contacts can be enabled/disabled with 
           Radou.params['save_contacts'].

           All proline residues and the N-terminal residue are skipped. See 
           calc_hbonds and calc_nh_contacts for optional kwargs.       

           Usage: PF()
       
           Returns: (array of residue indices, array of mean protection factors)"""
        # Setup residue/atom lists        
        hn_atms = Functions.extract_HN(self.t, log=self.params['logfile'])
        prolines = Functions.list_prolines(self.t, log=self.params['logfile'])
        # All protein residues except prolines
        reslist = [ r.index for r in self.top.residues if r.is_protein and r.index not in prolines[:,1] ]

        # Calc Nc/Nh
        hres, hbonds = self.calc_hbonds(hn_atms)
        cres, contacts = self.calc_nh_contacts(reslist)

        if not np.array_equal(hres, cres):
            raise Functions.HDX_Error("The residues analysed for Nc and Nh appear to be different. Check your inputs!")

        # Option to save outputs
        if self.params['save_contacts']:
            for i, residx in enumerate(hres):
                np.savetxt("Hbonds_%d.tmp" % self.top.residue(residx).resSeq, hbonds[i], fmt='%d') # Use residue indices internally, print out IDs
            for i, residx in enumerate(cres):
                np.savetxt("Contacts_%d.tmp" % self.top.residue(residx).resSeq, contacts[i], fmt='%d') # Use residue indices internally, print out IDs
        # Calc PF with phenomenological equation
        hbonds *= self.params['betah']     # Beta parameter 1
        contacts *= self.params['betac']   # Beta parameter 2
    
        pf = np.exp(hbonds + contacts)
        pf = np.mean(pf, axis=1)
        rids = np.asarray([ self.top.residue(i).resSeq for i in hres ])
        # Save PFs to separate log file, appending filenames for trajectories read as chunks
        if os.path.exists(self.params['outprefix']+"Protection_factors.dat"):
            filenum = len(glob.glob(self.params['outprefix']+"Protection_factors*"))
            np.savetxt(self.params['outprefix']+"Protection_factors_chunk_%d.dat" % (filenum+1), \
                       np.stack((rids, pf), axis=1), fmt=['%7d','%18.8f'], \
                       header="ResID  Protection factor") # Use residue indices internally, print out IDs
        else:    
            np.savetxt(self.params['outprefix']+"Protection_factors.dat", np.stack((rids, pf), axis=1), \
                       fmt=['%7d','%18.8f'], header="ResID  Protection factor") # Use residue indices internally, print out IDs

        return hres, pf

    def pro_omega_indices(self, prolines):
        """Calculates omega dihedrals (CA-C-N-CA) for all proline
           residues in a given prolines array from list_prolines.
    
           Usage: pro_omega_indices(prolines)
           Returns: (atom_indices, w_angles_by_frame)"""

        atom_names = ['CA', 'C', 'N', 'CA']
        offset = np.asarray([-1, -1, 0, 0])

        w_indices = np.zeros((len(prolines),4))
        for i, residx in enumerate(prolines[:,1]):
            curr_offset = offset + residx
            curr_atm_indices = []
            for x in zip(curr_offset, atom_names):
                # Cycle through previous CA/C, current N, CA
                curr_atm_indices.append(self.top.residue(x[0]).atom(x[1]).index)
            w_indices[i] = curr_atm_indices

        return w_indices, md.compute_dihedrals(self.t, w_indices)

    # Assignments for intrinsic rate adjustments:
    # 1) cis/trans prolines
    # 2) disulfides
    # 3) His protonation states
    # 4) N/C termini

    def assign_cis_proline(self):
        """Assigns cis-proline residues on a by-frame basis"""

        prolines = Functions.list_prolines(self.t, log=self.params['logfile'])
        outidxs, outangs = self.pro_omega_indices(prolines)
        for i, proidx in enumerate(prolines[:,1]):
            self.top.residue(proidx).cis_byframe = np.logical_and(outangs < np.pi/2, outangs > -1*np.pi/2)[:,i]
            if np.max(self.top.residue(proidx).cis_byframe) > 0:
                with open(self.params['logfile'], 'a') as f:
                    f.write("Cis-proline found at frame %d for residue %s!\n" % (np.argmax(self.top.residue(proidx).cis_byframe) + 1, self.top.residue(proidx).resSeq))

    def assign_disulfide(self):
        """Assigns residues involved in disulfide bridges"""

        # This selection syntax & assignment is untested
        sg = self.top.select('name SG and resname CYS or resname CYX')
        try:
            sg_coords = self.t.atom_slice(sg).xyz
        except IndexError: # Catch empty sg
            with open(self.params['logfile'], 'a') as f:
                f.write("No cysteines found in topology.\n")
            return
        self.top.create_standard_bonds()
        self.top.create_disulfide_bonds(sg_coords)
        # Assign disulfides (identical for each frame)
        for b in self.top._bonds:
            if all(i.element.symbol == 'S' for i in b):
                b[0].residue.disulf = np.ones(self.t.n_frames, dtype=bool)
                b[1].residue.disulf = np.ones(self.t.n_frames, dtype=bool)
                with open(self.params['logfile'], 'a') as f:
                    f.write("Disulfide found for residues %s - %s\n" \
                             % (b[0].residue, b[1].residue))

    def assign_his_protonation(self):
        """Assigns protonation state to HIS residues"""

        hisidx = [ r.index for r in self.top.residues if r.code == 'H' ]
        for i in hisidx:
            names = [ a.name for a in self.top.residue(i).atoms ]
            if all(n in names for n in ['HD1','HE2']): # Atom names for doubly protonated His (Hip)
                self.top.residue(i).HIP = np.ones(self.t.n_frames, dtype=bool)
                with open(self.params['logfile'], 'a') as f:
                    f.write("Protonated His assigned for residue %d\n" % self.top.residue(i).resSeq)
#            else:
#                self.top.residue(i).HIP = np.zeros(self.t.n_frames, dtype=bool)

    def assign_termini(self):
        """Assigns flags to N and C terminal residues"""

        for c in self.top.chains:
            c.residue(0).nterm = np.ones(self.t.n_frames, dtype=bool)
            c.residue(-1).cterm = np.ones(self.t.n_frames, dtype=bool)
            with open(self.params['logfile'], 'a') as f:
                f.write("N-terminus identified at: %s\nC-terminus identified at: %s\n" \
                         % (c.residue(0), c.residue(-1)))


    # Helper function to turn sequence-specific rate adjustments to intrinsic acid/base/water rates
    def _adj_to_rates(self, rate_adjs):
        """Helper function for kint().
           Calculates intrinsic rates for a given set of rate adjustments
           [ log(AL), log(BL), log(AR), log(BR) ] taken from Bai et al.

           Usage: _adj_to_rates(rate_adjs)
           Returns: intrinsic_rate"""

        # Calc reference rates at experimental temperature
        # / np.log(10) = conversion from ln to log10
        lgkAexp = self.params['kint_params']['lgkAref'] - (self.params['kint_params']['EaA'] \
                  /  np.log(10) / self.params['kint_params']['R']) * \
                  (1./self.params['kint_params']['Texp'] - 1./self.params['kint_params']['Tref'])
        lgkBexp = self.params['kint_params']['lgkBref'] - (self.params['kint_params']['EaB'] \
                  /  np.log(10) / self.params['kint_params']['R']) * \
                  (1./self.params['kint_params']['Texp'] - 1./self.params['kint_params']['Tref'])
        lgkWexp = self.params['kint_params']['lgkWref'] - (self.params['kint_params']['EaW'] \
                  /  np.log(10) / self.params['kint_params']['R']) * \
                  (1./self.params['kint_params']['Texp'] - 1./self.params['kint_params']['Tref'])

        # Calc log(kA||kB||kW)
        lgkA = lgkAexp + rate_adjs[0] + rate_adjs[2] - self.params['kint_params']['pD']
        lgkB = lgkBexp + rate_adjs[1] + rate_adjs[3] - self.params['kint_params']['pKD'] + self.params['kint_params']['pD']
        lgkW = lgkWexp + rate_adjs[1] + rate_adjs[3]

        kint = 10**lgkA + 10**lgkB + 10**lgkW
        #print(lgkAexp, lgkBexp, lgkWexp, 10**lgkA, 10**lgkB, 10**lgkW)
        return kint

# Intrinsic rate calc:
    def kint(self):
        """Function for calculating intrinsic rates of residues
           in a given topology
       
           Intrinsic exchange rates k_int are computed using equations below.
           k_int = k_A + k_B + k_W
           lgk_A = lgk_A,ref + lgA_L + lgA_R - pD
           lgk_B = lgk_B,ref + lgB_L + lgB_R - pOD
                 = lgk_B,ref + lgB_L + lgB_R - pK_D + pOD
           lgk_W = lgk_W,ref + lgB_L + lgB_R

           Default parameters for the above can be modified in the 
           Radou.params['kint_params'] dictionary. Sequence-based
           rate adjustments can be modified in the 'kint_adjs' and
           '_reordered_kint_adjs' dictionaries.

           Usage: kint()
           Returns: array of by-residue intrinsic rates"""

        kints = np.zeros(len(self.reslist))


     # Adjust residue names for: Cis-Pro, HIP, cystine bridges, GLH/ASH
        reslist = np.insert(self.reslist,0,self.reslist[0]-1) # Insert 'prev' residue for first index
        for i in reslist:
            curr = self.top.residue(i)
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
            if curr.name == 'GLU': # If Glu has a protonated carboxylate
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
        for chain in self.top.chains:
            chain.residue(0).name = 'NT'
            # Doesn't appead to be a standard name for COOH hydrogen - guess based on no. of bonds!
            try:
                if chain.residue(-1).atom('O').n_bonds > 1 or chain.residue(-1).atom('OXT').n_bonds > 1:
                    chain.residue(-1).name = 'CTH'
                    with open(self.params['logfile'], 'a') as f:
                        f.write("It looks like you have a neutral C-terminus (COOH) at residue %s?\n" % chain.residue(-1))
                else:
                    chain.residue(-1).name = 'CT'
            except KeyError:
                with open(self.params['logfile'], 'a') as f:
                    f.write("Residue %s at the end of a chain doesn't seem to have atoms named 'O'/'OXT'.\nI'm not treating it as a C-terminus.\n" % chain.residue(-1))

        reslist = np.delete(reslist, 0) # Remove i-1 residue we inserted above
        try:
            if np.array_equal(reslist, self.reslist):
                pass
            else:
                raise Functions.HDX_Error("Your residue lists for protection factors and intrinsic rates are different. Check your inputs!")
        except AttributeError:
            print("Please generate protection factors before running intrinsic rate calculations.")
            return
        for i, r in enumerate(reslist):
            curr = self.top.residue(r)
            prev = self.top.residue(r-1)
            # check for cispro
            if prev.name == 'PROC':
                with open(self.params['logfile'], 'a') as f:
                    f.write("Performing rate calculation on a by-frame basis for residue %s" % prev)
                adj_cis, adj_trans = self.params['_reordered_kint_adjs'][curr.name][0:2], self.params['_reordered_kint_adjs'][curr.name][0:2]
                adj_cis.extend(self.params['_reordered_kint_adjs']['PROC'][2:4])
                adj_trans.extend(self.params['_reordered_kint_adjs']['PRO'][2:4])
                kint_cis = self._adj_to_rates(adj_cis)
                kint_trans = self._adj_to_rates(adj_trans)
                kint_byframe = np.where(prev.cis_byframe, kint_cis, kint_trans)
                kints[i] = np.mean(kint_byframe) # Overall intrinsic rate is adjusted by mean population of cis-pro
            else:
                curr_adjs = self.params['_reordered_kint_adjs'][curr.name][0:2]
                curr_adjs.extend(self.params['_reordered_kint_adjs'][prev.name][2:4])
                kints[i] = self._adj_to_rates(curr_adjs)

            rids = np.asarray([ self.top.residue(i).resSeq for i in reslist ])
            # Save Kints to separate log file, appending filenames for trajectories read as chunks
        if os.path.exists(self.params['outprefix']+"Intrinsic_rates.dat"):
            filenum = len(glob.glob(self.params['outprefix']+"Intrinsic_rates*"))
            np.savetxt(self.params['outprefix']+"Intrinsic_rates_chunk_%d.dat" % (filenum+1), \
                       np.stack((rids, kints), axis=1), fmt=['%7d','%18.8f'], \
                       header="ResID  Intrinsic rate / min^-1 ") # Use residue indices internally, print out IDs
        else:    
            np.savetxt(self.params['outprefix']+"Intrinsic_rates.dat", np.stack((rids, kints), axis=1), \
                       fmt=['%7d','%18.8f'], header="ResID  Intrinsic rate / min^-1 ") # Use residue indices internally, print out IDs

        return kints

    # Deuterated fration by residue
    def dfrac(self, write=True):
        """Function for calculating by-residue deuterated fractions, for
           a set of Protection factors, intrinsic rates, and exposure times
           previously defined for the current Radou object.

           Usage: dfrac()
           Returns: [n_residues, n_times] 2D numpy array of deuterated fractions"""


        if len(set(map(len,[self.reslist, self.pfs, self.rates]))) != 1: # Check that all lengths are the same (set length = 1)
            raise Functions.HDX_Error("Can't calculate deuterated fractions, your residue/protection factor/rates arrays are not the same length.")

        try:
            fracs = np.zeros((len(self.reslist), len(self.params['times'])))
        except TypeError:
            fracs = np.zeros((len(self.reslist), 1))
        for i2, t in enumerate(self.params['times']):
            def _residue_fraction(pf, k, time=t):
                return 1 - np.exp(-k / pf * time)
            for i1, curr_frac in enumerate(itertools.imap(_residue_fraction, self.pfs, self.rates)):
                fracs[i1,i2] = curr_frac

        rids = np.asarray([ self.top.residue(i).resSeq for i in self.reslist ])
        # Write resfracs to separate file, appending filenames for trajectories read as chunks
        if write:
            if os.path.exists(self.params['outprefix']+"Residue_fractions.dat"):
                filenum = len(glob.glob(self.params['outprefix']+"Residue_fractions*"))
                np.savetxt(self.params['outprefix']+"Residue_fractions_chunk_%d.dat" % (filenum+1), \
                           np.hstack((np.reshape(rids, (len(self.reslist),1)), fracs)), \
                           fmt='%7d ' + '%8.5f '*len(self.params['times']), \
                           header="ResID  Deuterated fraction, Times / min: %s" \
                           % ' '.join([ str(t) for t in self.params['times'] ])) # Use residue indices internally, print out IDs
            else:    
                np.savetxt(self.params['outprefix']+"Residue_fractions.dat", \
                           np.hstack((np.reshape(rids, (len(self.reslist),1)), fracs)), \
                           fmt='%7d ' + '%8.5f '*len(self.params['times']), \
                           header="ResID  Deuterated fraction, Times / min: %s" \
                           % ' '.join([ str(t) for t in self.params['times'] ])) # Use residue indices internally, print out IDs
        return fracs

    def run(self, trajectory):
        """Runs a by-residue HDX prediction for the provided MDTraj trajectory

           Usage: run(traj)
           Returns: None (results are stored as Radou.resfracs)"""
        self.t = trajectory # Note this will add attributes to the original trajectory, not a copy
        self.n_frames = self.t.n_frames
        self.top = trajectory.topology.copy() # This does not add attributes to the original topology
        self.assign_cis_proline()
        self.assign_disulfide()
        self.assign_his_protonation()
        self.assign_termini()
        self.reslist, self.pfs = self.PF()
                                   
        self.rates = self.kint()
        self.resfracs = self.dfrac()
        print("Residue predictions complete")

### Add further classes for methods below here
