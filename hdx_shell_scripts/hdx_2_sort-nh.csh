#! /bin/csh
# Sort Nc-Nh results
# Author: Daniel Deredge <dderedge@rx.umaryland.edu>
# Adaptations by Anu Nagarajan/Richard Bradshaw

### N.B. This script is no longer necessary with gromacs 5 ###






###################################################################################
# When Nc-Nh are calculated, zero results are not output. Thus one should add     #
# missing outcomes before computing protection factors using the phenomenological #
# equation as the protection factors should be estimated as an ensemble average.  #
###################################################################################

###################
# Global settings #
###################

# Set variables
set dire = 1-nc-nh_G467                    # where you put output
set sres = 1                            # first residue ID (now no need to add 1)
set fres = 66                          # last  residue ID (now no need to add 1)
set ti   = 0                            # starting  time point of MD, ns
set tf   = 491.78                          # finishing time point of MD, ns
set dt   = 0.01                        # time interval, ns

# Make destination directory, if not exist
if ( ! -e $dire ) then
   mkdir $dire
endif
cd $dire

# Make an auxillary file
awk 'BEGIN{for(i='$ti';i<='$tf';i=i+'$dt') print i}' > tmp1 # print out time series
set nr = `wc -l tmp1 | awk '{print $1}'`                    # No. of timepoints (snapshots)

####################
# Sort Nc outcomes #
####################
if ( ! -e nc_old ) then
   mkdir nc_old
   mv nc_res*.dat nc_old
endif
rm -f nc_res*.dat
set i = $sres
while ( $i <= $fres )
    @ j = $i + 1 # proline resid's in MD
    set test = `grep -qx "$j" ../pro.list ; echo $status`
    if ( $test == 0 ) then
        goto nextnc # residue is proline, go to next residue
    else if ( $test == 1 ) then
        set k = `wc -l nc_old/nc_res$i.dat | awk '{print $1}'`
        if ( $k == $nr ) then
            cp nc_old/nc_res$i.dat .
            goto nextnc # there is no zero output, got to next residue
        else
            awk '{print $1}' nc_old/nc_res$i.dat > tmp2
            sort tmp1 tmp2 | uniq -u | awk '{print $0, 0}' > tmp3
            cp nc_old/nc_res$i.dat tmp4
            cat tmp3 >> tmp4
            sort -nk 1 tmp4 > nc_res$i.dat
            rm -f tmp[234]
        endif
    endif
    nextnc:
    @ i++
end

####################
# Sort Nh outcomes #
####################
if ( ! -e nh_old ) then
   mkdir nh_old
   mv nh_res*.dat nh_old
endif
rm -f nh_res*.dat
set i = $sres
while ( $i <= $fres )
    @ j = $i + 1 # proline resid's in MD
    set test = `grep -qx "$j" ../pro.list ; echo $status`
    if ( $test == 0 ) then
        goto nextnh # residue is proline, go to next residue
    else if ( $test == 1 ) then
        set k = `wc -l nh_old/nh_res$i.dat | awk '{print $1}'`
        if ( $k == $nr ) then
            cp nh_old/nh_res$i.dat .
            goto nextnh # there is no zero output, got to next residue
        else
            awk '{print $1}' nh_old/nh_res$i.dat > tmp2
            sort tmp1 tmp2 | uniq -u | awk '{print $0, 0}' > tmp3
            cp nh_old/nh_res$i.dat tmp4
            cat tmp3 >> tmp4
            sort -nk 1 tmp4 > nh_res$i.dat
            rm -f tmp[234]
        endif
    endif
    nextnh:
    @ i++
end

# Remove auxillary file
rm -f tmp1
