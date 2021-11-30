#!/bin/bash
cd BigBrainWarp
source scripts/init.sh

# economo atlas from fsaveage to BigBrain
dataDir=/data_/mica1/03_projects/casey/micasoft/parcellations/fsaverage7/
wd=/home/cpaq3/Desktop/DMN/utilities/
bigbrainwarp --in_space fsaverage --out_space bigbrain --wd $wd \
--in_lh $dataDir/lh.economo.annot --in_rh $dataDir/rh.economo.annot \
--desc economo --interp nearest

# Schaefer atlas from fsaveage to BigBrain
dataDir=/data_/mica1/03_projects/casey/micasoft/parcellations/fsaverage7/
wd=/home/cpaq3/Desktop/DMN/utilities/
bigbrainwarp --in_space fsaverage --out_space bigbrain --wd $wd \
--in_lh $dataDir/lh.Schaefer2018_400Parcels_7Networks_order.annot --in_rh $dataDir/rh.Schaefer2018_400Parcels_7Networks_order.annot \
--desc Schaefer2018_400Parcels_7Networks_order --interp nearest

# dmn probability map to fs_LR to fsaverage to BigBrain
dataDir=/data_/mica1/03_projects/casey/sandbox1/9_DMN/data/hcp/
bigbrainwarp --in_space fs_LR --out_space bigbrain --wd $dataDir \
    --in_lh $dataDir/tpl-fs_LR_hemi-L_desc-Kong_0.95.txt \
    --in_rh $dataDir/tpl-fs_LR_hemi-R_desc-Kong_0.95.txt \
    --desc Kong_0.95 --interp nearest

# data-driven cytoarchitectural axis of DMN from BigBrain to fsaverage
dataDir=/data_/mica1/03_projects/casey/sandbox1/9_DMN/output/
bigbrainwarp --in_space bigbrain --out_space fsaverage --wd $dataDir \
    --in_lh $dataDir/tpl-bigbrain_hemi-L_desc-DMN.txt \
    --in_rh $dataDir/tpl-bigbrain_hemi-R_desc-DMN.txt \
    --desc DMN --interp linear

# Kong parcellation to BigBrain
dataDir=/data_/mica1/03_projects/casey/sandbox1/9_DMN/data/hcp/
bigbrainwarp --in_space fs_LR --out_space bigbrain --wd $dataDir \
    --desc HCP_Kong_group --interp nearest \
    --in_lh $dataDir/tpl-fs_LR_hemi-L_desc-HCP_Kong_group.txt \
    --in_rh $dataDir/tpl-fs_LR_hemi-R_desc-HCP_Kong_group.txt
