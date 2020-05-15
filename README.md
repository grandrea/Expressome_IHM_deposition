This folder contains the cluster center model of the integrative modeling solution.



modeling.py : the modeling protocol employed with IMP 2.12

Folders:

cluster.0 : 
contains the .mrc files for the localisation probability densities for each region of interest
input: 
topology.txt: the integrative modeling set-up for IMP.
rnap_nusa.fasta: the sequences in the system
the restraints 
	density: RNAPNusA_job153_postprocess_masked.mrc 
	density GMM representation: gmm_700_Nov19.mrc
	DSS crosslinks: cl_all_nusg_rpoe_DSSO.csv
	DSSO crosslinks: cl_all_nusg_rpoe_DSS.csv

the starting models (*.pdb)
Sequence alignments for P75090 (delta) and P75591 (NusA). Other sequence alignments are found in the pdb file headers.

input/original_models:
the models coming out of comparative modeling, prior to coarse graining of flexible regions.


