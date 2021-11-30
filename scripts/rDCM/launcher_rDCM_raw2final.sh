#!/bin/bash

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 
# This bash script comprises the full analysis pipeline associated with the following paper:
#
#
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# -----------------------------------------
# Environment
# -----------------------------------------


# switch to new software stack (mainly to enable ghostscript)
. /cluster/apps/local/env2lmod.sh


# Load the ghostscript and matlab module
module load ghostscript/9.21
module load matlab


# output the new software stack (should have ghostscript and Matlab included)
module list


# create logs folder
mkdir logs



# -----------------------------------------
# Empirical analysis: whole-brain network
# -----------------------------------------


# Fit rs-fMRI data using rDCM for HCP dataset
for ID in {1..110}
do
	bsub -W 4:00 -J "job_rDCM_whole_HCP_$ID" matlab -singleCompThread -r "estimate_rDCM_fixed_MainExperiment(1,$ID,1,1)"
done


# Fit rs-fMRI data using rDCM for MICS dataset
for ID in {1..49}
do
	bsub -W 4:00 -J "job_rDCM_whole_MICS_$ID" matlab -singleCompThread -r "estimate_rDCM_fixed_MainExperiment(2,$ID,1,1)"
done


# Collect all results from the rDCM inversion
bsub -W 4:00 -J "job_get_rDCM_whole_HCP" -w 'ended("job_rDCM_whole_HCP*")' -oo "logs/log_rDCM_collect_whole_HCP" matlab -singleCompThread -r "get_rDCM_parameter_estimates_fixed_MainExperiment(1,1,1)"
bsub -W 4:00 -J "job_get_rDCM_whole_MICS" -w 'ended("job_rDCM_whole_MICS*")' -oo "logs/log_rDCM_collect_whole_MICS" matlab -singleCompThread -r "get_rDCM_parameter_estimates_fixed_MainExperiment(2,1,1)"


# Summarize the results from the whole-brain rDCM analysis
bsub -W 4:00 -J "job_summary_rDCM_whole_HCP" -w 'ended("job_get_rDCM_whole_HCP")' -oo "logs/log_rDCM_summary_whole_HCP" matlab -singleCompThread -r "summarize_rDCM_fixed_MainExperiment(1,1,1,0)"
bsub -W 4:00 -J "job_summary_rDCM_whole_MICS" -w 'ended("job_get_rDCM_whole_MICS")' -oo "logs/log_rDCM_summary_whole_MICS" matlab -singleCompThread -r "summarize_rDCM_fixed_MainExperiment(2,1,1,0)"


# Summarize the hierarchy results from the whole-brain rDCM analysis
bsub -W 4:00 -J "job_hierarchy_rDCM_whole_HCP" -w 'ended("job_summary_rDCM_whole_HCP")' -oo "logs/log_rDCM_hierarchy_whole_HCP" matlab -singleCompThread -r "hierarchy_rDCM_fixed_MainExperiment(1,1,1)"
bsub -W 4:00 -J "job_hierarchy_rDCM_whole_MICS" -w 'ended("job_summary_rDCM_whole_MICS")' -oo "logs/log_rDCM_hierarchy_whole_MICS" matlab -singleCompThread -r "hierarchy_rDCM_fixed_MainExperiment(2,1,1)"
