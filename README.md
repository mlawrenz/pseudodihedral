pseudodihedral
==============

example code to compute pseudo dihedrals for groups of atoms in a simulation. 

In this simple test case I used a single frame pdb (new-aln-state847.cent.pdb), and generate 4 groups composed of 2 residues (all atoms). I compute the COM for each of the 4 groups, and pass it into the msmbuilder dihedral calculation code (see the import from msmbuilder) to get the dihedral.

I checked the value I computed with VMD by doing the following:
* Load in the pdb, set four selections with the 2 residues, and get the COM coordinate:
set sel [ atomselectop top "resid 1 2"] 
puts measure center $sel weight mass
 {6.23124 3.43234 4.35556}

* Copy these four COM coordinates into test.pdb
I just replaced the coordinates of 4 atoms (which are meaningless here) with the COM coordinates from VMD, and made sure the format is all correct so VMD can load it. Then I use the VMD measure dihedral command to double check my result.


