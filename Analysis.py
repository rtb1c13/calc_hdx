#!/usr/bin/env python

# Analysis/plotting functions for HDX analysis

import Functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as correl

class Analyze():
    """Class to contain analysis & plotting methods for HDX predictions"""

    def __init__(self, resobj, top, **extra_params):
        """Initialises Analyze object from a Method object with by-residue results"""
        try:
            self.reslist = resobj.reslist
            self.resfracs = resobj.resfracs
            self.n_frames = resobj.t.n_frames
            self.top = top
        except AttributeError:
            raise Functions.HDX_Error("Error when copying results from prediction to analysis objects - have you made any HDX predictions yet?")
        self.params = extra_params



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
           and the segment/average deuteration information to "Segment_average_fractions.tmp"

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

                idxs = np.where(np.logical_and( self.reslist > start, self.reslist <= end ))[0] # > s[0] skips the first residue in segment
                aves[i1, i2] = np.mean(self.resfracs[idxs, i2])

        np.savetxt("Segment_average_fractions.tmp", np.hstack((self.segres, aves)), \
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
            self.correls[idx] = correl(self.expfracs[:,idx], self.segfracs[:,idx])
            self.MSE[idx] = np.mean(self.expfracs[:,idx] - self.segfracs[:,idx])
            self.MUE[idx] = np.mean(np.abs(self.expfracs[:,idx] - self.segfracs[:,idx]))
            
    def run(self, figs=False):
        """Runs a by-segment HDX prediction and optionally graphs results"""

        self.read_segfile()
        self.segfracs = self.segments(self.top)
        if self.params['expfile'] is not None:
            self.read_expfile()
            self.desc_stats()
            


### Add further classes below here








