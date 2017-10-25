#! /bin/csh
# Compute deuterium fraction vs. time for a segment generated in HDX experiment
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# Do NOT count prolines, see Radou_BiophysJ_2014_v107_p983                        #
###################################################################################
# Before using the script, make sure you have                                     #
# pro.list - listing prolines                                                     #
# exp_seg.list - listing first & last residue ID's of segments in HDX experiment  #
###################################################################################
# The script can NOT treat terminal residues, need to update in the future.       #
###################################################################################

###################
# Global settings #
###################

# Set variables
set dire = "7-d-fraction-seg" # where you put output

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif

##############################
# Compute Deuterium fraction #
##############################
cd $dire
rm -f df_seg_*.dat # remove old files
set nseg = `wc -l ../exp_seg.list | awk '{print $1}'` # No. of segments in HDX experiment
@ seg = 1
while ( $seg <= $nseg )
    set sres = `awk '{if(NR=='$seg') print $1}' ../exp_seg.list` # first resid of the segment
    set fres = `awk '{if(NR=='$seg') print $2}' ../exp_seg.list` # last  resid of the segment
    @ sres = $sres - 1 # Segres index shifts from resids used in MD
    @ fres = $fres - 1 # Segres index shifts from resids used in MD
    @ i = $sres  # Current segres index
    while ( $i <= $fres )
        @ j = $i + 1 # proline checks go back fo resids, not idxs
        set test = `grep -qx "$j" ../pro.list ; echo $status` # check if res$i is PRO
        if ( $test == 0 ) then # in case that res$i is PRO, got to next residue
            echo "Skipping residue "$i
            goto nextres
        else if ( $test == 1 ) then # res$i is not PRO
            if ( $i == $sres ) then # This misses the case where a starting residue is a PRO
                awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ../6-d-fraction/df_res$i.dat > tmp
            else
                awk '{print $2}' ../6-d-fraction/df_res$i.dat | awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' >> tmp
            endif
        endif
        nextres:
        @ i++
    end
    @ sres = $sres + 1 # Segres index returns to resids as per MD
    @ fres = $fres + 1 # Segres index returns to resids as per MD
    awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' tmp | gawk '{sum=0;for(i=2;i<=NF;i++){sum+=$i};sum/=(NF-1);printf "%4.5f %1.6f \n", $1, sum}' > df_seg_${sres}-$fres.dat
    rm -f tmp
    @ seg++
end
