#!/usr/bin/env python

# Analysis/plotting functions for HDX analysis

import Functions
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from scipy.stats import pearsonr as correl

class Analyze():
    """Class to contain analysis & plotting methods for HDX predictions"""

    def __init__(self, resobj, top, **extra_params):
        """Initialises Analyze object from a Method object with by-residue results"""
        try:
            self.reslist = resobj.reslist
            self.resfracs = resobj.resfracs
            self.pfs = resobj.pfs
            self.n_frames = resobj.t.n_frames
            self.top = top
        except AttributeError:
            raise Functions.HDX_Error("Error when copying results from prediction to analysis objects - have you made any HDX predictions yet?")
        self.params = resobj.params
        try:
            self.params.update(extra_params)
        except (TypeError, ValueError):
            print("Couldn't load extra parameters for analysis (maybe they weren't provided?).\nUsing previous parameters from %s object." % resobj)
                   

    def read_segfile(self):

        # segfile should contain at most 2 columns, startres & endres
        self.segres = np.loadtxt(self.params['segfile'], dtype='i8')  # ResIDs will be converted to indices with dictionary in segments function

    def read_expfile(self):
        """Reads an experimental data file for comparison to predicted data.

           Experimental results file should be formatted as follows:
           Seg_start   Seg_end   Time_1   Time_2   Time_3 ... [Time_n]

           This is the same format as the printout of predicted results"""

        # Check I'm not loading in too many timepoints
        try:
            expt = np.loadtxt(self.params['expfile'], dtype=[ ('segres', np.int32, (2,)),\
                              ('fracs', np.float64, (len(self.params['times']),)) ])
        except ValueError:
            raise Functions.HDX_Error("""There's a problem with your experimental data file, it's shorter than the number of timepoints you evaluated""")
        # Now check I'm not loading in too few timepoints
        try:
            expt = np.loadtxt(self.params['expfile'], dtype=[ ('segres', np.int32, (2,)),\
                              ('fracs', np.float64, (len(self.params['times']) + 1,)) ])
            raise Functions.HDX_Error("""There's a problem with your experimental data file, it's longer than the number of timepoints you evaluated""")
        except ValueError:
            pass
        # Check expt = predicted
        if np.array_equal(self.segres, expt['segres']):
            self.expfracs = expt['fracs']
        else:
            raise Functions.HDX_Error("The experimental segments read from %s and predicted segments read from %s don't match!" % (self.params['segfile'], self.params['expfile']))
                                       
    
    def segments(self, top):
        """Function to average residue deuterated fractions over
           a given set of peptide segments. The first residue in each
           segment will not be included in the averaging, as it is assumed
           to be 100% back-exchanged during analysis.
    
           Residue indices provided in the given list are converted to 
           residue IDs from the given trajectory's topology. Currently this
           remumbering will only work for single chain topologies with sequential
           numbering. If a residue in a segment is not found (e.g. a truncated
           N/C terminus), the next residue is chosen as the start/end point instead.
 
           Writes info on skipped residues to logfile "HDX_analysis.log" by default
           and the segment/average deuteration information to "Segment_average_fractions.dat"

           Usage: segments(traj, reslist, fracs, segfile_name, times, [ log="HDX_analysis.log" ])
           Returns: [n_segs, 2] 2D numpy array of segment start/end residue IDs, 
                [n_segs, n_times] 2D numpy array of segment deuterated fractions at each timepoint"""

        res2idx = {}
        with open(self.params['logfile'], 'a') as f:
            f.write("Now converting residue numbers to indices for segment averaging:\n")
        for idx, res in enumerate(top.residues):
            if res.is_protein:
#                res2idx[str(res.chain.index) + '-' + str(res.resSeq)] = idx
                res2idx[res.resSeq] = idx # Only works for single chain or sequential numbers, no re-use of resnums
            else:
                with open(self.params['logfile'], 'a') as f:
                    f.write("Skipping residue: %s, not a protein residue\n" % res)
        
        
        self.read_segfile()
        try:
            aves = np.zeros((len(self.segres), len(self.params['times'])))
        except TypeError:
            aves = np.zeros((len(self.segres), 1))
            self.params['times'] = [self.params['times']]    

        # Info for 'skip_first'
        if self.params['skip_first']:
            for i1, s in enumerate(self.segres):
                with open(self.params['logfile'], 'a') as f:
                    f.write("'Skip_first' is set. Not including residue %s in averaging for segment %s-%s.\n" \
                            % (top.residue(res2idx[s[0]]), s[0], s[1]))
        else:
            for i1, s in enumerate(self.segres):
                with open(self.params['logfile'], 'a') as f:
                    f.write("'Skip_first' is NOT set. Including residue %s in averaging for segment %s-%s.\n" \
                            % (top.residue(res2idx[s[0]]), s[0], s[1]))

        for i2, t in enumerate(self.params['times']):
            for i1, s in enumerate(self.segres):
                try:
                    start = res2idx[s[0]]
                except KeyError:
                    with open(self.params['logfile'], 'a') as f:
                        f.write("Didn't find residue %s in protein. Using residue %s as startpoint instead.\n" \
                             % (s[0], top.residue(0)))
                    start = 0
                try:
                    end = res2idx[s[1]]
                except KeyError:
                    with open(self.params['logfile'], 'a') as f:
                        f.write("Didn't find residue %s in protein. Using residue %s as endpoint instead.\n" \
                                 % (s[0], top.residue(-1)))
                    end = top.residue(-1).index

                if self.params['skip_first']:
                    idxs = np.where(np.logical_and( self.reslist > start, self.reslist <= end ))[0] # > start skips
                else:
                    idxs = np.where(np.logical_and( self.reslist >= start, self.reslist <= end ))[0] # >= start incs
                    
                aves[i1, i2] = np.mean(self.resfracs[idxs, i2])

        np.savetxt(self.params['outprefix']+"Segment_average_fractions.dat", np.hstack((self.segres, aves)), \
                   fmt='%6d %6d ' + '%8.5f '*len(self.params['times']), header="Res 1   Res2  Times / min: %s" \
                   % ' '.join([ str(t) for t in self.params['times'] ]))
        with open(self.params['logfile'], 'a') as f:
            f.write("Segment averaging complete.\n")

        return aves


    def desc_stats(self):
        """Calculates descriptive statistics of segments compared to expt
           for all timepoints"""

        self.correls = np.zeros(len(self.params['times']))
        self.MUE = np.zeros(len(self.params['times']))
        self.MSE = np.zeros(len(self.params['times']))
        for idx, t in enumerate(self.params['times']):
            self.correls[idx] = correl(self.segfracs[:,idx], self.expfracs[:,idx])[0]
            self.MSE[idx] = np.mean(self.segfracs[:,idx] - self.expfracs[:,idx])
            self.MUE[idx] = np.mean(np.abs(self.segfracs[:,idx] - self.expfracs[:,idx]))
        with open(self.params['outprefix']+"Descriptive_statistics.dat", 'a') as f:
            np.savetxt(f, self.correls, header="Pearson's R correlation, Times / min: %s" \
                       % ' '.join([ str(t) for t in self.params['times'] ]), fmt='%8.6f')
            np.savetxt(f, self.MSE, header="Mean signed error / frac, Pred. - Expt., Times / min: %s" \
                       % ' '.join([ str(t) for t in self.params['times'] ]), fmt='%10.8f')
            np.savetxt(f, self.MUE, header="Mean unsigned error / frac, Pred. - Expt., Times / min: %s" \
                       % ' '.join([ str(t) for t in self.params['times'] ]), fmt='%10.8f')

    def print_summaries(self):
        """Print summary PF and resfrac results - for example of a method
           object that has had results summed over multiple chunks."""

        
        # Save PFs to 'SUMMARY' file
        try:
            rids = np.asarray([ self.top.residue(i).resSeq for i in self.reslist ])
            if os.path.exists(self.params['outprefix']+"SUMMARY_protection_factors.dat"):
                filenum = len(glob.glob(self.params['outprefix']+"SUMMARY_protection_factors*"))
                np.savetxt(self.params['outprefix']+"SUMMARY_protection_factors_%d.dat" % (filenum+1), \
                           np.stack((rids, self.pfs), axis=1), fmt=['%7d','%18.8f'], \
                           header="ResID  Protection factor") # Use residue indices internally, print out IDs
            else:    
                np.savetxt(self.params['outprefix']+"SUMMARY_protection_factors.dat", np.stack((rids, self.pfs), axis=1), \
                           fmt=['%7d','%18.8f'], header="ResID  Protection factor") # Use residue indices internally, print out IDs
        except AttributeError:
            print("Can't write summary protection factors - perhaps you haven't calculated them yet?")
            return
        try:
            if os.path.exists(self.params['outprefix']+"SUMMARY_residue_fractions.dat"):
                filenum = len(glob.glob(self.params['outprefix']+"SUMMARY_residue_fractions*"))
                np.savetxt(self.params['outprefix']+"SUMMARY_residue_fractions_%d.dat" % (filenum+1), \
                           np.hstack((np.reshape(rids, (len(self.reslist),1)), self.resfracs)), \
                           fmt='%7d ' + '%8.5f '*len(self.params['times']), \
                           header="ResID  Deuterated fraction, Times / min: %s" \
                           % ' '.join([ str(t) for t in self.params['times'] ])) # Use residue indices internally, print out IDs
            else:    
                np.savetxt(self.params['outprefix']+"SUMMARY_residue_fractions.dat", \
                           np.hstack((np.reshape(rids, (len(self.reslist),1)), self.resfracs)), \
                           fmt='%7d ' + '%8.5f '*len(self.params['times']), \
                           header="ResID  Deuterated fraction, Times / min: %s" \
                           % ' '.join([ str(t) for t in self.params['times'] ])) # Use residue indices internally, print out IDs
        except AttributeError:
            print("Can't write summary residue fractions - perhaps you haven't calculated them yet?")
            return
        
    
    def run(self, figs=False):
        """Runs a by-segment HDX prediction and optionally graphs results"""

        self.read_segfile()
        self.segfracs = self.segments(self.top)
        if self.params['expfile'] is not None:
            self.read_expfile()
            self.desc_stats()
            


### Add further classes below here








