#! /bin/csh
# Compute deuterium fraction vs. time for each residue
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# Deuterium fraction time series of a residue is calculated use the equation      #
# D(t)=1-exp[(-k_int/P-factor)*t]                                                 #
# See Radou_BiophysJ_2014_v107_p983 for details                                   #
###################################################################################
# Before using the script, make sure you have                                     #
# pro.list - listing prolines                                                     #
###################################################################################
# The script can NOT treat terminal residues, need to update in the future.       #
###################################################################################

###################
# Global settings #
###################

# Set variables
set dire = 6-d-fraction                 # where you put output
set sres = 1                            # first residue ID
set fres = 512                            # penultimate residue ID (last isn't calculated as C-terminal)

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif

##############################
# Compute Deuterium fraction #
##############################
cd $dire
rm -f * # remove old files
@ i = $sres # start from initial res index (2nd residue)
while ( $i < $fres )
    @ j = $i + 1 # proline resids in MD
    set test = `grep -qx "$j" ../pro.list ; echo $status` # check if res$i is PRO
    if ( $test == 0 ) then # in case that res$i is PRO, got to next residue
        goto nextres
    else if ( $test == 1 ) then # res$i is not PRO
        # get mean P-factor, remembering pf has redidue idxs, not IDs
        set pf = `grep -v "#" ../4-p-factor/pf_res$i.dat | awk '{sum+=$4}END{printf "%12.5f", sum/NR}'`
        # get pertinent intrinsic exchange rate, remembering k_iexr has residue IDs, not idxs
        set kint = `awk '{if($1=='$i'+1) printf "%10.10f", $3}' ../5-kint/k_iexr.dat`
        echo $i $kint $pf >> df_plot.dat
        foreach t ( 0 0.167 1 10 120 ) # time in min, MUST match experimental values!
            set df = `awk 'BEGIN{printf "%1.6f", 1-exp(-'$kint'/'$pf'*'$t')}'` # deuterium fraction
            echo $t $df | awk '{printf "%4.5f %1.6f \n", $1, $2}' >> df_res$i.dat
        end
    endif
    nextres:
    @ i++
end
