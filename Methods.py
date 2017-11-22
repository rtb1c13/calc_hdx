#!/usr/bin/env python

# Class for HDX trajectories, inherited from MDTraj

import mdtraj as md
import numpy as np

class Radou():
    """Class for Radou-style analysis"""

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
                        'outprefix' : None }
        self.params.update(extra_params) # Update main parameter set from kwargs

        # Default rate adjustments for adjacent amino acidsi
        # Each key is a residue name, each value is [ lgAL, lgAR, lgBL, lgBR ]
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

    def calc_contacts(self, qidx, cidx, cutoff):
        """Calculates contacts between 'query' and 'contact' atom selections
           within a specified cutoff (default = 0.65, for coordinates in nm).
           Periodicity is included in MDtraj function by default.
           Usage: calc_contacts(qidx, cidx, cutoff).

           Qidx and cidx are the atom index lists to search for contacts from
           and to respectively (e.g. from amide NH to all heavy atoms).

           Returns count of contacts for each frame in supplied trajectory."""

        try:
            byframe_ctacts = md.compute_neighbors(self.t, cutoff, qidx, haystack_indices=cidx)
        except TypeError:
            print("Now calculating contacts to single atom, idx %d" % qidx)
            qidx = np.array([qidx])
            byframe_ctacts = md.compute_neighbors(self.t, cutoff, qidx, haystack_indices=cidx)
        return map(lambda x: len(x), byframe_ctacts)

    def _calc_hbonds_contacts(self, HN):
        """Calculates number of protein H-bonds for a particular atom index
           using the 'contacts' method. Bonds to all protein O* or N* evaluated
           by default, optionally all non-protein too (including waters) if 
           Radou.params['protonly'] is True.
       
           Usage: _calc_hbonds_contacts(atom)"""

        # Get N index in same residue as current HN atom
        getN4H = lambda _: self.top.atom(_).residue.atom('N').index
        if self.params['protonly']:
            c = self.top.select("protein and (symbol O or symbol N) and not index %s" % getN4H(HN))
        else:
            c = self.top.select("all (symbol O or symbol N) and not index %s" % getN4H(HN))

        hbond_counts = self.calc_contacts(HN, c, self.params['cut_Nh'])
        return hbond_counts

    def _calc_hbonds_bh(self, HN, minfreq=0.0, **kwargs):
        """Calculates number of protein H-bonds for a particular atom index
           using the 'Baker-Hubbard' method. Default donor-acceptor distance < 0.25 nm
           + angle > 120 degrees.
           Reports all H-bonds (minimum freq=0.0) by default. Bonds to all protein 
           O* or N* evaluated by default, optionally all non-protein too 
           (including waters) if Radou.params['protonly'] is True. For other kwargs
           see mdtraj.geometry.hbond._get_bond_triplets or ..._compute_bounded_geometry.
       
           Usage: _calc_hbonds_bh(atom, [minfreq, **kwargs])
           Returns: n_frames length array of H-bond counts for desired atom"""

        # Atoms for H-bonds includes protein or all O*, N* and single HN hydrogen

        if self.params['protonly']:
            c = self.t.atom_slice(self.top.select("protein and (symbol O or symbol N) or index %s" % HN))
        else:
            c = self.t.atom_slice(self.top.select("all (symbol O or symbol N) or index %s" % HN))

        # Call internal functions of md.baker_hubbard directly to return
        # distances & angles, otherwise only bond_triplets averaged across
        # a trajectory are returned
        bond_triplets = md.geometry.hbond._get_bond_triplets(c.topology, **kwargs)
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

    def calc_hbonds(self, donors, **kwargs):
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

        if self.params['skip_first']:
            donors = donors[1:] # Remove first atom/residue from list

        try:
            total_counts = np.zeros((len(donors), self.t.n_frames))
        except TypeError:
            total_counts = np.zeros((1, self.t.n_frames))
        for i, v in enumerate(donors):
            total_counts[i] = methods[self.params['hbond_method']](v, **kwargs)

        reslist = [ self.top.atom(a).residue.index for a in donors ]
#        hbonds = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), total_counts), axis=1) # Array of [[ Res idx, Contact count ]]

        return np.asarray(reslist), total_counts

    def calc_nh_contacts(reslist):
        """Calculates contacts between each NH atom and the surrounding heavy atoms,
           excluding those in residues n-2 to n+2.
    
           By default contacts < 0.65 nm are calculated, and only protein-heavys,
           are included, but can include all heavys if desired. Also skips first 
           residue (N-terminus) in a residue list by default too.

           Usage: calc_nh_contacts(traj, reslist, [cutoff=0.65, skip_first=True, protein=True])
           Returns: n_res x n_frames 2D array of contacts per frame for each residue in reslist"""

        # Check if current atom is a heavy atom
        is_heavy = lambda _: self.top.atom(_).element.symbol is not 'H'

        if self.params['skip_first']:
            reslist.pop(0) # Remove first atom/residue from list

        contact_count = np.zeros((len(reslist), self.t.n_frames))
        for idx, res in enumerate(reslist):
            robj = self.top.residue(res)
            excl_idxs = range(robj.index - 2, robj.index + 3, 1) # Exclude n-2 to n+2 residues

### *** ### Calls external (select_residxs)
            inv_atms = select_residxs(self.t, excl_idxs, protonly=self.params['protonly'], invert=True) # At this stage includes H + heavys
### *** ###
            heavys = inv_atms[ np.array( [ is_heavy(i) for i in inv_atms ] ) ] # Filter out non-heavys

            contact_count[idx] = self.calc_contacts(robj.atom('N').index, heavys, cutoff=self.params['cut_Nc'])

#        contacts = np.concatenate((np.asarray([reslist]).reshape(len(reslist),1), contact_count), axis=1) # Array of [[ Res idx, Contact count ]]
        return np.asarray(reslist), contact_count

    def PF(self, **kwargs):
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
### *** ### calls extract_HN/list_prolines
        # Setup residue/atom lists        
        hn_atms = extract_HN(self.t)
        prolines = list_prolines(self.t)
        # All protein residues except prolines
        reslist = [ r.index for r in self.top.residues if r.is_protein and r.index not in prolines[:,1] ]

        # Calc Nc/Nh
        hres, hbonds = self.calc_hbonds(hn_atms, **kwargs)
        cres, contacts = self.calc_nh_contacts(reslist)

        if not np.array_equal(hres, cres):
            raise HDX_Error("The residues analysed for Nc and Nh appear to be different. Check your inputs!")

        # Option to save outputs
        if self.params['save_contacts']:
            for i, residx in enumerate(hres):
                np.savetxt("Hbonds_%d.tmp" % self.top.residue(residx).resSeq, hbonds[i], fmt='%d') # Use residue indices internally, print out IDs
            for i, residx in enumerate(cres):
                np.savetxt("Contacts_%d.tmp" % self.top.residue(residx).resSeq, contacts[i], fmt='%d') # Use residue indices internally, print out IDs
        # Calc PF with phenomenological equation
        hbonds *= self.params['betah']        # Beta parameter 1
        contacts *= self.params['betac']   # Beta parameter 2
    
        pf = np.exp(hbonds + contacts)
        pf = np.mean(pf, axis=1)
        rids = np.asarray([ self.top.residue(i).resSeq for i in hres ])
        # Save PFs to separate log file
        np.savetxt(self.params['outprefix']+"Protection_factors.tmp", np.stack((rids, pf), axis=1), \
                   fmt=['%7d','%18.8f'], header="ResID  Protection factor") # Use residue indices internally, print out IDs

        return hres, pf




    def run(trajectory):
        """Runs a by-residue HDX prediction for the provided trajectory"""
        self.t = trajectory # Note this will add attributes to the original trajectory, not a copy
        self.top = trajectory.topology.copy() # This does not add attributes to the original topology
        assign_cis_proline(self.t, log=self.params['logfile'])
        assign_disulfide(self.t, log=self.params['logfile'])
        assign_his_protonation(self.t, log=self.params['logfile'])
        assign_termini(self.t, log=self.params['logfile'])
        self.reslist, self.pfs = PF(self.t, self.params['hbond_method'], self.params['save_contacts'],)
                                   
        k, r = kint(protcurr, hres)
        fracs += dfrac(protlipcurr, hres, pf, k, times)
