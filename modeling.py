"""
#############################################
##  IMP Tutorial Script
##
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: Papers where this data is shown...
#
"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import IMP.bayesianem
import IMP.bayesianem.restraint


import os
import sys

# fix for mpi error
#import DLFCN as dl
#sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)


import IMP.mpi

#---------------------------
# Define Input Files
#---------------------------
datadirectory = "/scratch/graziade/imp_segmented_map_feb20/input/"
topology_file = datadirectory + "/topology.txt"
target_gmm_file = datadirectory + '/gmm_700_Nov19.txt'

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames_XL = 30000
if '--test' in sys.argv: num_frames_XL = 100

num_mc_steps = 10
if '--test' in sys.argv: num_mc_steps = 15

#--------------------------
# Create movers
#--------------------------

# rigid body movement params
rb_max_trans = 6.00
rb_max_rot = 0.1

# flexible bead movement
bead_max_trans = 4.00

#floppy bodies max trans
# not used
fb_max_trans = 5.0

################################################
#

#--------------------------------
# Build the Model Representation
#--------------------------------

# Initialize model
m = IMP.Model()


# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file, fasta_dir=datadirectory,
                                           pdb_dir=datadirectory, gmm_dir=datadirectory)


bm = IMP.pmi.macros.BuildSystem(m, force_create_gmm_files=True)
bm.add_state(topology)

root_hierarchy, dof = bm.execute_macro(max_rb_trans=rb_max_trans, 
        max_rb_rot=rb_max_rot, 
        max_bead_trans=bead_max_trans, 
        max_srb_rot=rb_max_rot,
        max_srb_trans=rb_max_trans)

outputobjects = []
# sampleobjects = [] # not used anywhere

#outputobjects.append(hierarchy)
# sampleobjects.append(root_hierarchy)



fixed_particles = []
for prot in ["Q50301", 
        "P75560", 
        "P75581", 
        "P41205", 
        "P46775", 
        "Q50304", 
        "P75179", 
        "P75049",
        "longRNAR1",
        "longRNAR2",
        "DNA1N", 
        "DNA1T",
        "P75271",
        "30SsubunitE",
        "30SsubunitF",
        "30SsubunitJ",
        "30SsubunitK",
        "30SsubunitL",
        "30SsubunitM",
        "30SsubunitN",
        "30SsubunitO",
        "30SsubunitP",
        "30SsubunitR",
        "30SsubunitT",
        "30Ssubunitj"]:
    fixed_particles += IMP.atom.Selection(root_hierarchy, molecule=prot).get_selected_particles()

##block NTD of alpha
mysel=IMP.atom.Selection(root_hierarchy, molecule="Q50295.1", residue_indexes=range(0,250))
fixed_particles += mysel.get_selected_particles()
mysel2=IMP.atom.Selection(root_hierarchy, molecule="Q50295.2", residue_indexes=range(0,250))
fixed_particles += mysel2.get_selected_particles()
#
#block core of  beta
mysel3=IMP.atom.Selection(root_hierarchy, molecule="P78013", residue_indexes=range(0,984))
fixed_particles += mysel3.get_selected_particles()
mysel4=IMP.atom.Selection(root_hierarchy, molecule="P78013", residue_indexes=range(1005,1392))
fixed_particles += mysel4.get_selected_particles()
#
#mysel7=IMP.atom.Selection(root_hierarchy, molecule="P75271", residue_indexes=range(614,1270))
#fixed_particles += mysel7.get_selected_particles()

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads, fixed_rbs = dof.disable_movers(fixed_particles,
                                            [IMP.core.RigidBodyMover,
                                             IMP.pmi.TransformMover])

IMP.pmi.tools.shuffle_configuration(root_hierarchy,
                                     max_translation=400,
                                     excluded_rigid_bodies=fixed_rbs,
                                     verbose=True,
                                     cutoff=8.0,
                                     niterations=500)

# IMP.pmi.tools.shuffle_configuration(root_hierarchy)

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files


# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hierarchy, 
                                         resolution=10,
                                         kappa=2.0)
ev.add_to_model()


# ----------- tutorial ------------
# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hierarchy)
for mol in mols:
    molname = mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)
# -----------------

outputobjects.append(ev)

dof.optimize_flexible_beads(200)

# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)

from IMP.pmi.io.crosslink import FilterOperator as FO
cldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("fromSite")
cldbkc.set_residue2_key("ToSite")


print("DSS BRO")
cldb1 = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb1.create_set_from_file(datadirectory + "cl_all_nusg_rpoe_DSS.csv")

cl_dss = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hierarchy,
                                                                            CrossLinkDataBase=cldb1,
                                                                            length=25.0,
                                                                            slope=0.02,
                                                                            resolution=1.0,
                                                                            filelabel='DSS',
                                                                            label="DSS",
                                                                            weight=5)

cl_dss.set_sigma_is_sampled(True)
cl_dss.set_psi_is_sampled(True)
cl_dss.add_to_model()             # crosslink must be added to the model
# sampleobjects.append(xl1) #crosslink restraint is storing a sampled particle
outputobjects.append(cl_dss)
dof.get_nuisances_from_restraint(cl_dss)
print("DSSO, DUDE")

cldb2 = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb2.create_set_from_file(datadirectory + "cl_all_nusg_rpoe_DSSO.csv")

cl_dsso = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hierarchy,
                                                                            CrossLinkDataBase=cldb2,
                                                                            length=25.0,
                                                                            slope=0.02,
                                                                            resolution=1.0,
                                                                            filelabel='DSSO',
                                                                            label="DSSO",
                                                                            weight=10)

cl_dsso.set_sigma_is_sampled(True)
cl_dsso.set_psi_is_sampled(True)
cl_dsso.add_to_model()             # crosslink must be added to the model
outputobjects.append(cl_dsso)
dof.get_nuisances_from_restraint(cl_dsso)


#print("HMW maaaaan")

#cldb3 = IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
#cldb3.create_set_from_file(datadirectory + "HMWfrac_cl_all_nusg_wo_rb.csv")
#
#cl_hmwfrac = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hierarchy,
#                                                                            CrossLinkDataBase=cldb3,
#                                                                            length=30.0,
#                                                                            slope=0.02,
#                                                                            resolution=1.0,
#                                                                            filelabel='HMWfrac',
#                                                                            label="HMWfrac",
#                                                                            weight=1)
#
#cl_hmwfrac.set_sigma_is_sampled(False)
#cl_hmwfrac.add_to_model()             # crosslink must be added to the model
#outputobjects.append(cl_hmwfrac)
#dof.optimize_flexible_beads(100)
####

# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-Bayesian

sel = IMP.atom.Selection(hierarchy=root_hierarchy, representation_type=IMP.atom.DENSITIES)
densities=sel.get_selected_particles()

gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                 target_fn=target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=0.01,
                                                 target_radii_scale=3.0,
                                                 target_is_rigid_body=False)
gem.add_to_model()
gem.set_label("EM_Score")
outputobjects.append(gem)





#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hierarchy,
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    output_objects=outputobjects,
                                    crosslink_restraints=[cl_dss, 
                                        cl_dsso],    # allows XLs to be drawn in the RMF files
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=False,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=20.0,
                                    simulated_annealing_minimum_temperature_nframes=1000,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=20.0,
                                    number_of_best_scoring_models=0,
                                    replica_exchange_swap=True,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames_XL,
                                    test_mode=False,
                                    save_coordinates_mode="25th_score")

# Start Sampling
mc1.execute_macro()

#print("gaussian restraint")
#print(gem.evaluate())
print("cl_dss restraint")
print(cl_dss.evaluate())
#print("cl_hmw restraint")
#print(cl_hmwfrac.evaluate())
print("cl_dsso restraint")
print(cl_dsso.evaluate())
print("ev restraint")
print(ev.evaluate())

