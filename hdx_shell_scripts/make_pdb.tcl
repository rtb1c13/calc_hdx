# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)
#
# Read in psf of full protein
mol new {/u/lucy/anu/Leut/from_pacific/3TT3/MD/P419/y268a/3tt3_y268a_full.psf} type {psf} first 0 last -1 step 1 waitfor 1
# Read in crd (Charmm coords) of full system, energy minimised
mol addfile {/u/lucy/anu/Leut/from_pacific/3TT3/MD/P419/y268a/3tt3_y268a_full.crd} type {cor} first 0 last -1 step 1 waitfor 1 0
# Write PDB of full system
animate write pdb {/u/bradshaw/WORK/LeuT/Anu_files/3TT3_E290_K288_Y168/Run_3/Leut_y268a_min.pdb} beg 0 end 0 skip 1 0
# Write PDB of protein only
animate write pdb {/u/bradshaw/WORK/LeuT/Anu_files/3TT3_E290_K288_Y168/Run_3/Leut_desolv.pdb} beg 0 end 0 skip 1 sel [atomselect top protein]
#
quit
# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)
# end of log file.
