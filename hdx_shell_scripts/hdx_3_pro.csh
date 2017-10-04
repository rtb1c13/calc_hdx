#! /bin/csh -x
# Determine peptide bond conformation for prolines
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# As the intrinsic exchange rates when proline neighbors the residues differs     #
#    among cis- & trans-PRO, one needs to determine the peptide bond conformation #
#    i.e., compute omega dihedrals for prolines.                                  #
# Check Bai_Proteins_1993_v17_p75 for details                                     #

###################################################################################

module load gromacs/5.1.4

###################
# Global settings #
###################

# Set variables
set dire = 3-pro_omega                # where you put output
set tpr  = ../Leut_HDX.tpr            # run input file
set trj  = Leut_stripped.trr          # trajectory file

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif

########################################
# Compute omega dihedrals for prolines #
########################################
rm -f $dire/*
foreach i (`cat pro.list`)
    @ j = $i - 1 # Prev residue for w angle
    echo "r $j & a CA" > tmp.txt
    echo "r $j & a C" >> tmp.txt
    echo "r $i & a N" >> tmp.txt
    echo "r $i & a CA" >> tmp.txt
    echo "10 | 11 | 12 | 13" >> tmp.txt
    echo "name 14 PRO_res${j}_omega" >> tmp.txt
    echo "del 10-13" >> tmp.txt
    echo "q" >> tmp.txt
    gmx make_ndx -f $tpr -o res$j.ndx < tmp.txt
    echo "10" > tmp.txt
    gmx angle -f $trj -n res$j.ndx -type dihedral -ov res$j.xvg < tmp.txt
    sed '/#/d' res$j.xvg | sed '/@/d' > $dire/res$j.dat
    rm -f tmp.txt res$j.* angdist.xvg
    set cisp = `awk '{if(($2<=90)&&($2>=-90)) print $2}' $dire/res$j.dat | wc -l`
    if ( $cisp > 0 ) then
        echo "PRO$j has cis conformation!!!!!" # in case of cis-proline
    endif
end
