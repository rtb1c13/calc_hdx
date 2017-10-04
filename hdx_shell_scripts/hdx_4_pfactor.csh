#! /bin/csh
# Compute protection factors of residues
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# To correlate simulation data with experimental results, the phenomenological    #
# equation (Vendruscolo_JAmChemSoc_2003_v125_p15686) is used:                     #
# ln(P-factor)=parc*Nc+parh*Nh                                                    #
# Two variables should be computed.                                               #
# Nc: No. of contacts between heavy atoms of OTHER residues to the amide hydrogen #
# Nh: No. of hydrogen bonds to the amide hydrogen                                 #
# Parameters parc & parh were optimized in Best_Structure_2006_v14_p97            #
###################################################################################

###################
# Global settings #
###################

# Set variables
set dire = 4-p-factor-hbonds            # where you put output
set parc = 0.35                         # beta_c in the phenomenological equation
set parh = 2                            # beta_h in the phenomenological equation
set sres = 1                            # first residue ID
set fres = 512                          # last  residue ID

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif

#############################
# Compute Protection factor #
#############################
cd $dire
rm -f pf_res*.dat # remove old file
set i = $sres
while ( $i <= $fres )
    @ j = $i + 1 # proline resids in MD
    set test = `grep -qx "$j" ../pro.list ; echo $status`
    if ( $test == 0 ) then
        goto nextres # residue is proline, go to next residue
    else if ( $test == 1 ) then
        echo "# Time Nc Nh Pfactor" > pf_res$i.dat
        cp ../1-nc-nhbonds/nc_res$i.dat tmp1
        awk '{print $2}' ../1-nc-nhbonds/nh_res$i.dat > tmp2
        paste tmp1 tmp2 | awk '{printf "%4.4f %3d %3d %12.5f \n", $1, $2, $3, exp('$parc'*$2+'$parh'*$3)}' >> pf_res$i.dat
        rm -f tmp*
    endif
    nextres:
    @ i++
end
