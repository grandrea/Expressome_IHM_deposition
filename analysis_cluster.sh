#!/bin/bash

#SBATCH -J IMP_analysis	# Job Name
#SBATCH --nodes=1 		# Anzahl Knoten N
#SBATCH --ntasks-per-node=20 	
#SBATCH --mem=32G              # 500MiB resident memory pro node

##Max Walltime vorgeben:
#SBATCH --time=18:00:00 # Erwartete Laufzeit

#Auf Standard-Knoten rechnen:
#SBATCH --partition=standard

#Job-Status per Mail:
#SBATCH --mail-type=NONE
#SBATCH --mail-user=andrea.graziadei@tu-berlin.de


#conda init bash
source /home/users/g/graziadei/anaconda3/etc/profile.d/conda.sh
source /home/users/g/graziadei/.bashrc


conda activate imp-realease

cd /scratch/graziadei/imp_segmented_map_feb20/scripts_sampcon5_ranges

#rm cl_headers.txt
#
#python ../imp-sampcon/pyext/src/select_good_scoring_models.py -rd ../ -rp "run" -sl "CrossLinkingMassSpectrometryRestraint_Distance_" -pl  CrossLinkingMassSpectrometryRestraint_Data_Score_DSSO CrossLinkingMassSpectrometryRestraint_Data_Score_DSS ExcludedVolumeSphere_None GaussianEMRestraint_EM_Score Total_Score ConnectivityRestraint_P75271 ConnectivityRestraint_P75591 -alt 0.9 -aut 1.0 -mlt 0.0 -mut 30.0
#
#python ../imp-sampcon_Feb20/pyext/src/select_good_scoring_models.py -rd ../ -rp run -sl "CrossLinkingMassSpectrometryRestraint_Distance_" "GaussianEMRestraint_EM_Score" "ExcludedVolumeSphere_None" \
#        -pl CrossLinkingMassSpectrometryRestraint_Data_Score_DSSO CrossLinkingMassSpectrometryRestraint_Data_Score_DSS ExcludedVolumeSphere_None GaussianEMRestraint_EM_Score Total_Score ConnectivityRestraint_P75271 ConnectivityRestraint_P75591 \
#	-alt 0.90 0.0 0.0 \
#	-aut 1.0 5775.0 60.0 \
#	-mlt 0.0 0.0 0.0  \
#	-mut 35.0 0.0 0.0 -e
#
#cp -rv ../good_scoring_models .

#python ../imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py -n rnap_nusa_ribo -p ../good_scoring_models/ -g 1.0 -d density_ranges_all.txt -m cpu_omp -c 4 -dt 10.0

#python ../imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py -n rnap_nusa_ribo -p ./good_scoring_models/ -g 0.5 -d density_ranges_all.txt -m cpu_omp -c 24 -dt 20.0

#python ../imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py -n rnap_nusa_ribo -p ../good_scoring_models/ -g 0.1 -d density_ranges_all.txt -m cpu_omp -c 24 --subunit "P75591"
#mkdir all_proteins
#cp -rv * all_proteins
#
python ../imp-sampcon_Feb20_all/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py -n rnap_nusa_ribo -p ./good_scoring_models -g 0.5 -d density_ranges_all.txt -m cpu_omp -c 20 --subunit "P75591"

