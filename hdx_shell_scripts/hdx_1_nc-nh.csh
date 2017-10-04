#! /bin/csh
# Compute Nh and Nc for computing deuterium fraction
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

###################################################################################
# To correlate simulation data with experimental results, the phenomenological    #
# equation (Vendruscolo_JAmChemSoc_2003_v125_p15686) is used:                     #
# ln(P-factor)=parc*Nc+parh*Nh                                                    #
# Two variables should be computed.                                               #
# Nc: No. of contacts between heavy atoms of OTHER residues to the amide hydrogen #
# Nh: No. of hydrogen bonds to the amide hydrogen                                 #
# Nc cutoff: 6.5 Angstrom, amide nitrogen and heavy atoms of OTHER residues       #
# Nh cutoff: 2.4 Angstrom, amide hydrogen and acceptor oxegen/nitrogen            #
# Cutoffs from Best_Structure_2006_v14_p97 & Radou_BiophysJ_2014_v107_p983        #
###################################################################################

module load gromacs/5.1.4

###################
# Global settings #
###################

# Set variables
set trj  = Leut_stripped.trr            # trajectory file #CHANGE
set tpr  = ../Leut_HDX.tpr              # run input file #CHANGE
set pdb  = ../Leut_desolv.pdb           # structure file #CHANGE
set dire = 1-nc-nh                        # where you put output
set sres = 2                            # first residue ID (1st residue is NT-ACE in GROMACS)
set fres = 513                          # last  residue ID (1st residue is NT-ACE in GROMACS)
set cutc = 0.65                         # cutoff for Nc, nm
set cuth = 0.24                         # cutoff for Nh, nm

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir -p $dire
endif

# List prolines ##### CHANGE THIS PART TO PICK THE RIGHT CATEGORY : ORIGINAL had 10; Changed to 14 as that was the right one for my system #ANU
gmx make_ndx -f $pdb -o tmp.ndx << EOF
r PRO
splitres 10
del 0-10
q
EOF

grep 'PRO' tmp.ndx | awk '{print $2}' | sed 's/PRO_PRO_//g' > pro.list
rm -f tmp.ndx

# List residues for later analysis
awk 'NR>1 {print $6, $4}' $pdb | uniq > res.list # Assumes first line is CRYST record
#sed -i 's/HSD/HIS/g' res.list
sed -i 's/HSE/HSD/g' res.list # convert for script 5


##############
# Compute Nc #
##############
rm -f $dire/nc_res*.dat # remove old files
set i = $sres
while ( $i <= $fres )
    set test = `grep -qx "$i" pro.list ; echo $status`
    if ( $test == 0 ) then
        goto nextnc # residue is proline, go to next residue
    else if ( $test == 1 ) then
        @ pa = $i - 1
        @ pb = $i - 2
#        @ pc = $i - 3
        @ aa = $i + 1
        @ ab = $i + 2
#        @ ac = $i + 3
        # Make *.ndx for res$i
        echo "t NH1 & a N & r $i" > tmp.txt                    # amide nitrogen of res$i
        echo "name 10 amide_n_res$i" >> tmp.txt		### CHANGE the number for gromacs index based on your system
#        echo "! r 1 $i & 2" >> tmp.txt                         # heavy atoms of all except residue i
#        echo "! r 1 $pa $i $aa & 2" >> tmp.txt                 # heavy atoms of all except i, i+/-1
        echo "! r 1 $pb $pa $i $aa $ab & 2" >> tmp.txt         # heavy atoms of all except i, i+/-1/2
#        echo "! r 1 $pc $pb $pa $i $aa $ab $ac & 2" >> tmp.txt # heavy atoms of all except i, i+/-1/2/3
        echo "name 11 heavy_no-neighbor" >> tmp.txt	#### CHANGE the number for gromacs index based on your system
        echo "q" >> tmp.txt
        gmx make_ndx -f $tpr -o res$i.ndx < tmp.txt  # make *.ndx for res$i. MUST use *.tpr here
#        echo "10" > tmp.txt
#        echo "11" >> tmp.txt
        # Calcualte Nc for res$i
        gmx select -f $trj -s $tpr -n res$i.ndx -select 'group 11 and within '$cutc' of group 10' -oi tmp.dat
        @ j = $i - 1 # get real resid
        # NR/100 = time in ns for 10 ps intervals
        awk '{print NR/100, $2}' tmp.dat > $dire/nc_res$j.dat
        rm -f tmp.* res$i.ndx
    endif
    nextnc:
    @ i++
end

##############
# Compute Nh #
##############
rm -f $dire/nh_res*.dat # remove old files
set i = $sres
while ( $i <= $fres )
    set test = `grep -qx "$i" pro.list ; echo $status`
    if ( $test == 0 ) then
        goto nextnh # residue is proline, go to next residue
    else if ( $test == 1 ) then
        # Make *.ndx for res$i
        echo "t H & a HN & r $i" > tmp.txt        # define amide hydrogen of res$i
        echo "name 10 amide_h_res$i" >> tmp.txt
        echo "t NH1 & a N & r $i" >> tmp.txt      # define amide nitrogen of res$i
        echo "t N* O* & ! 11" >> tmp.txt          # define O/N of all residues except amide N of res$i
        echo "del 11" >> tmp.txt
        echo "name 11 nitrogen_oxygen" >> tmp.txt
        echo "q" >> tmp.txt
        gmx make_ndx -f $tpr -o res$i.ndx < tmp.txt  # make *.ndx for residue $i. MUST use *.tpr here
#        echo "10" > tmp.txt
#        echo "11" >> tmp.txt
        # Calcualte Nc for res$i
        gmx select -f $trj -s $tpr -n res$i.ndx -select 'group 11 and within '$cuth' of group 10' -oi tmp.dat 
        @ j = $i - 1 # get real resid
        # NR/100 = time in ns for 10 ps intervals
        awk '{print NR/100, $2}' tmp.dat > $dire/nh_res$j.dat
        rm -f tmp.* res$i.ndx
    endif
    nextnh:
    @ i++
end

# Move pro.list to hdx for future use
#mv pro.list hdx/
