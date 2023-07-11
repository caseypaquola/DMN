#!/bin/bash
cd BigBrainWarp
source scripts/init.sh

# data-driven cytoarchitectural axis of DMN from BigBrain to fsaverage
dataDir=~/GitHub/DMN/output/
bigbrainwarp --in_space bigbrain --out_space fsaverage --wd $dataDir \
    --in_lh $dataDir/tpl-bigbrain_hemi-L_desc-DMN.txt \
    --in_rh $dataDir/tpl-bigbrain_hemi-R_desc-DMN.txt \
    --desc DMN --interp linear