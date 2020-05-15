
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
# import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.pmi.mmcif
import ihm

import os
import sys

# fix for mpi error
#import DLFCN as dl
#sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)


import IMP.mpi

#---------------------------
# Define Input Files
#---------------------------
datadirectory = "input/"
#datadirectory = "/home/andrea/Dropbox/imp_runs/imp_rnap_nusA_Sept19_Nterm_placed_longRNA_NusG/input/"
topology_file = datadirectory + "/topology_swismo.txt"
target_gmm_file = datadirectory + '/gmm_700_Nov19.txt'
analysis_directory = './cluster.0/'

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
if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(None)
    bm.system.add_protocol_output(po)
    po.system.title = "M. pneumoniae expressome"
    s = po.system
    s.software[0].version = '2.12.0'
    s.software[1].version = '2.12.0'
    s.software.append(ihm.Software(name='SWISS-MODEL', classification='comparative modeling',
                                   description='a fully automated protein structure homology-modelling server',
                                   location='https://swissmodel.expasy.org/', version='2019-11-21'))
    s.authors = ['OReilly FJ', 'Xue L', 'Graziadei A', 'Sinn L',
                 'Lenz S', 'Tegunov D', 'Bloetz C', 'Hagen WJH', 'Cramer P', 'Stuelke J', 'Mahamid J', 'Rappsilber J'
                 ]
    # Add publication
#    po.system.citations.append(ihm.Citation.from_pubmed_id(25161197))

bm.dry_run = '--dry-run' in sys.argv

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

####


# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-Bayesian
#em_components = IMP.pmi.tools.get_densities(root_hierarchy)
#gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components, target_gmm_file, scale_target_to_mass=True, slope=0.000001, weight=180.0)
#gemt.add_to_model()
#outputobjects.append(gemt)

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
                                    test_mode=bm.dry_run,
                                    )

# Start Sampling
mc1.execute_macro()

last_step = s.orphan_protocols[-1].steps[-1]
last_step.num_models_end = 14400000 # TODO calculate properly

if '--mmcif' in sys.argv:
    # import ihm.cross_linkers
    import ihm.location
    # import ihm.model
    import RMF
    import IMP.rmf

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 100 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=14400000,
                            num_models_end=20765))

    mg = ihm.model.ModelGroup(name="Cluster 0")
    # Add to last state
    po.system.state_groups[-1][-1].append(mg)
    e = ihm.model.Ensemble(model_group=mg,
                           num_models=20131,
                           post_process=analysis.steps[-1],
                           name="Cluster 0")
    po.system.ensembles.append(e)

    # # # Create an ensemble for the cluster (warning: _add_simple_ensemble
    # # # is subject to change in future IMP releases) and deposit a single
    # # # representative model (let's say it's frame 42 from the output RMF file)
    # e = po._add_simple_ensemble(analysis.steps[-1],
    #                             name="Cluster 0", num_models=1,
    #                             drmsd=1, num_models_deposited=1,
    #                             localization_densities={}, ensemble_file=None)
    # e = po.set_ensemble_file()
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only(analysis_directory + 'cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [root_hierarchy])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    # del rh
    model = po.add_model(e.model_group)


    template_info_dict = {'Q50295.0': {'PDB': '6FLQ', 'PDB_chain': 'A', 'PDB_seq_range': (1, 329), 'seq_range': (29, 325), 'seq_id': 0.41},
                          'Q50295.1': {'PDB': '6FLQ', 'PDB_chain': 'B', 'PDB_seq_range': (1, 329), 'seq_range': (29, 325), 'seq_id': 0.41},
                          'P78013.0': {'PDB': '6FLQ', 'PDB_chain': 'C', 'PDB_seq_range': (1, 1342), 'seq_range': (14, 1356), 'seq_id': 0.41},
                          'P75271.0': {'PDB': '6FLQ', 'PDB_chain': 'D', 'PDB_seq_range': (1, 1407), 'seq_range': (14, 1280), 'seq_id': 0.41},
                          'P75049.0': {'PDB': '6C6U', 'PDB_chain': 'N', 'PDB_seq_range': (1, 181), 'seq_range': (13, 146), 'seq_id': 0.22},
                          'Q50301.0': {'PDB': '3J9W', 'PDB_chain': 'AE', 'PDB_seq_range': (1, 166), 'seq_range': (63, 219), 'seq_id': 0.516},
                          'P75560.0': {'PDB': '3J9W', 'PDB_chain': 'AB', 'PDB_seq_range': (1, 246), 'seq_range': (20, 244), 'seq_id': 0.353},
                          'P75581.0': {'PDB': '3J9W', 'PDB_chain': 'AJ', 'PDB_seq_range': (1, 102), 'seq_range': (7, 107), 'seq_id': 0.406},
                          'P41205.0': {'PDB': '3J9W', 'PDB_chain': 'AC', 'PDB_seq_range': (1, 218), 'seq_range': (2, 215), 'seq_id': 0.431},
                          'P46775.0': {'PDB': '3J9W', 'PDB_chain': 'AD', 'PDB_seq_range': (1, 200), 'seq_range': (2, 202), 'seq_id': 0.457},
                          'Q50304.0': {'PDB': '3J9W', 'PDB_chain': 'AH', 'PDB_seq_range': (1, 132), 'seq_range': (10, 142), 'seq_id': 0.508},
                          'P75179.0': {'PDB': '3J9W', 'PDB_chain': 'AI', 'PDB_seq_range': (1, 130), 'seq_range': (1, 132), 'seq_id': 0.569},

                          }
    starting_models = [segment.starting_model for segment in model.representation]
    for start_m_i in starting_models:
        if start_m_i is not None:
            unit_name = start_m_i.asym_unit.asym.details
            if unit_name in ['P75591.0', 'P75090.0', 'DNA1N.0', 'DNA1T.0', 'longRNAR2.0', 'longRNAR1.0'] or '30Ssubunit' in unit_name:
                continue
            start_m_i.templates = [ihm.startmodel.Template(
                dataset=ihm.dataset.PDBDataset(ihm.location.PDBLocation(template_info_dict[unit_name]['PDB'])),
                                    asym_id=template_info_dict[unit_name]['PDB_chain'], # asym_id (chain) of template (pdb)
                                    seq_id_range=template_info_dict[unit_name]['PDB_seq_range'],
                                    template_seq_id_range=template_info_dict[unit_name]['seq_range'],
                                    sequence_identity=template_info_dict[unit_name]['seq_id'])]
             #
    # add localisation densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    # analysis_directory += 'LPD_'
    asym = po.asym_units['P75049.0']
    # Add path to a local output file
    loc = ihm.location.OutputFileLocation(analysis_directory + 'LPD_P75049_loop.mrc')
    den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
    # Add to ensemble
    e.densities.append(den)

    asym = po.asym_units['P75090.0']
    # Add path to a local output file
    loc = ihm.location.OutputFileLocation(analysis_directory + 'LPD_P75090.mrc')
    den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
    # Add to ensemble
    e.densities.append(den)

    asym = po.asym_units['Q50295.0']
    # Add path to a local output file
    loc = ihm.location.OutputFileLocation(analysis_directory + 'LPD_Q50295_CTD.mrc')
    den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
    # Add to ensemble
    e.densities.append(den)

    asym = po.asym_units['P75591.0']
    # Add path to a local output file
    for domain_mrc in ['LPD_P75591_CTD.mrc', 'LPD_P75591_KH.mrc', 'LPD_P75591_S1.mrc', 'LPD_P75591_N.mrc']:
        loc = ihm.location.OutputFileLocation(analysis_directory + domain_mrc)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
        # Add to ensemble
        e.densities.append(den)

    asym = po.asym_units['P75271.0']
    # Add path to a local output file
    for domain_mrc in ['LPD_P75271_loop1.mrc', 'LPD_P75271_loop2.mrc']:
        loc = ihm.location.OutputFileLocation(analysis_directory + domain_mrc)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
        # Add to ensemble
        e.densities.append(den)

    asym = po.asym_units['P78013.0']
    # Add path to a local output file
    for domain_mrc in ['LPD_P78013_loop1.mrc', 'LPD_P78013_loop2.mrc']:
        loc = ihm.location.OutputFileLocation(analysis_directory + domain_mrc)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
        # Add to ensemble
        e.densities.append(den)

# em, = [r for r in s.restraints
#        if isinstance(r, ihm.restraint.EM3DRestraint)]
# d = em.dataset
# d.parents[0].location = ihm.location.EMDBLocation('EMD-10680')

for r in s.restraints:
    if isinstance(r, ihm.restraint.CrossLinkRestraint):
        if 'DSSO' in r.dataset.location.path:
            r.dataset.location = ihm.location.PRIDELocation('PXD017711')
        elif 'DSS' in r.dataset.location.path:
            r.dataset.location = ihm.location.PRIDELocation('PXD017695')

# repo = ihm.location.Repository(doi="xxx", root="./", url="https://github.com/Rappsilber-Laboratory/M.-pneumoniae-expressome/")
# s.update_locations_in_repositories([repo]) #


if '--mmcif' in sys.argv:
    po.finalize()
    with open('expressome_pdbdev.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])
#     po.flush()

    os.system('cp expressome_pdbdev.cif expressome_pdbdev.ihm')
