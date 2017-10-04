#!/bin/csh
# Compute intrinsic exchange rates for each residue in a protein
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# Intrinsic exchange rates k_int are computed using equations below.              #
# k_int = k_A + k_B + k_W                                                         #
# lgk_A = lgk_A,ref + lgA_L + lgA_R - pD                                          #
# lgk_B = lgk_B,ref + lgB_L + lgB_R - pOD                                         #
#       = lgk_B,ref + lgB_L + lgB_R - pK_D + pOD                                  #
# lgk_W = lgk_W,ref + lgB_L + lgB_R                                               #
# lgk_A,ref/lgk_B,ref/lgk_W,ref: check Table I & III of Bai_Proteins_1993_v17_p75 #
#                    T-dependence described by Van 't Hoff equation               #
# lgA_L/lgA_R/lgB_L/lgB_R: check Table II of Bai_Proteins_1993_v17_p75            #
# pK_D: T-dependent, check Covington_JPhysChem_1966_v70_p3820                     #
# Check Bai_Proteins_1993_v17_p75 for details.                                    #
###################################################################################
# Before using the script, make sure you have following three lists               #
# res.list - listing resnames & resids; HSE/HSD ---> HIS, leave HSP as they are.  #
#            If pKa's for His unknown, use HIS for all.                           #
#            Modify protonated GLU to GLUP                                        #
# pro.list - listing prolines                                                     #
# hdx_sd-exrate_new.dat - copied from Table II of Bai_Proteins_1993_v17_p75       #
###################################################################################
# Now the scripts is incapable of treating cis-PRO, but it does NOT cause error.  #
# The script can NOT treat terminal residues, need to update in the future.       #
###################################################################################

###################
# Global settings #
###################

# Set variables
set dire = 5-kint-PDLAref               # where you put output
set sres = 1                            # first residue ID
set fres = 513                          # last  residue ID
set Texp = 298                          # experimental temperature, K
set pKD  = 14.87                        # pKa of deuterium water at Texp
set pD   = 7.4                          # experimental pD

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif

################################################
# Compute k_A, k_B, k_W at your experimental T #
################################################
set lgkAref = 1.62     # lg reference reaction rates of  acid-catalyzed HDX, /M/min (low salt)
set lgkBref = 10.05    # lg reference reaction rates of  base-catalyzed HDX, /M/min (low salt)
set lgkWref = -1.5     # lg reference reaction rates of water-catalyzed HDX, /M/min (low salt)
set EaA     = 14       # activation energy for  acid-catalyzed HDX, kcal/mol
set EaB     = 17       # activation energy for  base-catalyzed HDX, kcal/mol
set EaW     = 19       # activation energy for water-catalyzed HDX, kcal/mol
set R       = 0.001987 # gas constant, kcal/mol/K
set Tref    = 293      # reference temperature, K

set lgkAexp = `echo "scale=5; $lgkAref-$EaA/2.303/$R*(1/$Texp-1/$Tref)" | bc`
set lgkBexp = `echo "scale=5; $lgkBref-$EaB/2.303/$R*(1/$Texp-1/$Tref)" | bc`
set lgkWexp = `echo "scale=5; $lgkWref-$EaW/2.303/$R*(1/$Texp-1/$Tref)" | bc`

echo $lgkAexp $lgkBexp $lgkWexp

####################################
# Compute intrinsic exchange rates #
####################################
cd $dire
rm -f k_iexr.dat # remove old file
# Work on non-terminal residues first
@ i = $sres + 1 # start from 2nd residue
while ( $i < $fres )
    @ j = $i  # Proline resids (not idxs), same as in MD
    set test = `grep -qx "$j" ../pro.list ; echo $status` # check if res$i is PRO
    if ( $test == 0 ) then
        goto nextres # res$i is PRO, go to next residue
    else if ( $test == 1 ) then # res$i is not PRO
        @ j = $i - 1 # resid of the residue to the left  of res$i
        @ k = $i + 1 # resid of the residue to the right of res$i
        # get resnames for residues j-i-k
        set resni = `awk '{if($1=='$i') print $2}' ../res.list` # resname of res$i
        set resnj = `awk '{if($1=='$j') print $2}' ../res.list` # resname of res$j
        set resnk = `awk '{if($1=='$k') print $2}' ../res.list` # resname of res$k
        echo $resnj $resni $resnk
        # determine lgAL, lgBL
        ### Modified by RTB, 15-Sep-17, to calculate NH to the left of sidechain $i, not $j
        ### See Bai et al. 1993
###        if ( $resnk == "PRO" ) then # res$k is PRO; res$i is to the left of PRO
###            set lgAL = 0.00
###            set lgBL = 0.00
###        else if ( $resnk == "HSD" ) then # res$k is HSD/HSE; res$i is to the left of HSE/HSD
###            set lgAL = 0.00
###            set lgBL = `grep HIS ../hdx_sd-exrate_new.dat | awk '{print $4}'`
        if ( $resni == "HSD" ) then # res$i is HSD/HSE; NH is to the left of HSE/HSD
            set lgAL = 0.00 # Is this correct? Strictly this is N/D, not 0.00
            set lgBL = `grep HIS ../hdx_sd-exrate_new.dat | awk '{print $4}'`
        else # res$i is not PRO/HSE/HSD
            set lgAL = `awk '{if($1=="'$resni'") print $2}' ../hdx_sd-exrate_new.dat`
            set lgBL = `awk '{if($1=="'$resni'") print $4}' ../hdx_sd-exrate_new.dat`
        endif
        # determine lgAR, lgBR
        if ( $resnj == "PRO" ) then # res$j is PRO; res$i is to the right of PRO
            set lgAR = `grep PRO ../hdx_sd-exrate_new.dat | sed '2d' | awk '{print $3}'` # trans-PRO
            set lgBR = `grep PRO ../hdx_sd-exrate_new.dat | sed '2d' | awk '{print $5}'` # trans-PRO
        else if ( $resnj == "HSD" ) then # res$j is HSD/HSE; res$i is to the right of HSE/HSD
            set lgAR = 0.00 # Again, strictly this is N/D, not 0.00
            set lgBR = `grep HIS ../hdx_sd-exrate_new.dat | awk '{print $5}'`
        else # res$j is not PRO/HSE/HSD
            set lgAR = `awk '{if($1=="'$resnj'") print $3}' ../hdx_sd-exrate_new.dat`
            set lgBR = `awk '{if($1=="'$resnj'") print $5}' ../hdx_sd-exrate_new.dat`
        endif
        # compute lgk_A/lgk_B/lgk_W
        set lgkA = `echo "scale=10; $lgkAexp+$lgAL+$lgAR-$pD" | bc`
        set lgkB = `echo "scale=10; $lgkBexp+$lgBL+$lgBR-$pKD+$pD" | bc`
        set lgkW = `echo "scale=10; $lgkWexp+$lgBL+$lgBR" | bc`
        # compute k_int
        set kA = `awk -v a=$lgkA 'BEGIN{printf "%3.10f", exp(2.303*a)}'`
        set kB = `awk -v a=$lgkB 'BEGIN{printf "%3.10f", exp(2.303*a)}'`
        set kW = `awk -v a=$lgkW 'BEGIN{printf "%3.10f", exp(2.303*a)}'`
        set kint = `echo "scale=10; $kA+$kB+$kW" | bc`
        echo $lgkA $lgkB $lgkW $kA $kB $kW $kint
        echo $i $resni $kint | awk '{printf "%4s %3s %10.10f \n", $1, $2, $3}' >> k_iexr.dat
    endif
    nextres:
    @ i++
end
# Now deal with terminal residues
