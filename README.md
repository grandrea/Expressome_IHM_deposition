This folder contains the cluster center model of the integrative modeling solution.

Files:

expressome.ihm :

The integrative model calculated by IMP. Can be opened in chimerax and is be deposited in PDB-Dev.
Once opened, contains all information on crosslinks and density used to generate the model, as well as the actual model and the localisation probability densities of each domain(in "Results").
The model is calculated based on DSS and DSSO crosslinks, and the multibody-refined RNAP-NusA density (found in "input").

Crosslinks from the RNAP alpha subunit (Q50295) appear four times, as the protein is present in two copies. Nevertheless, violations are computed on the shortest link.
N.B.: The "starting models" contained in the ihm are not aligned in the correct coordinate frame of the solution. As such, they do not represent the integrative modeling result. The position of each domain in the result is found in the "Results" tab or can be viewed by opening Modeling_cluster_center_NoCTD1.pdb .

Modeling_cluster_center.pdb:
The pdb version of the ihm file. 
-All coarse grained regions are removed. 


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


