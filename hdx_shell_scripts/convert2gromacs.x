#!/bin/bash

# Converts input files to gromacs format and makes tpr using standard mdp file
# Makes use of mdconvert utility from MDTraj

module load gromacs/5.1.4

# Filenames etc.
INTRJ=/u/lucy/anu/Leut/from_pacific/leut-f108y/total.dcd
OUTTRJ=HDX_test/Leut_stripped.trr
LASTATOMIDX=8257   # Zero indexed

# 1) Create atom index and strip DCD
seq -s " " 0 $LASTATOMIDX > atom_indices.dat
#mdconvert $INTRJ -o $OUTTRJ -a atom_indices.dat 

# 2) Create PDB (update make_pdb.tcl with filenames first)
vmd -dispdev text -e inputs/make_pdb.tcl

# Gromacs to create gro, tpr files (interactive)
gmx pdb2gmx -f Leut_desolv.pdb -o Leut_desolv.gro -p Leut_desolv.top -glu -ter 
# Uncomment if you need to check group numbers
#make_ndx -f Leut_desolv.gro -o Leut_desolv.ndx
gmx grompp -f inputs/dummy.mdp -c Leut_desolv.gro -p Leut_desolv.top -o Leut_HDX.tpr 
