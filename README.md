THE DEFAULT MODE NETWORK
----------------------------------------------

![Alt text](figures/figure_methods_overview.png?raw=true "Title")

Welcome to the code repository of "The Unique Cytoarchiecture and Wiring of the Human Default Mode Network". If you would like to dig deeper into the methods or reproduce aspects of the analyses, then please follow through these scripts.


### 👩‍💻️ Key scripts
**cytoarchitecture.m**
- Evaluates the heterogeneity of the DMN according to cortical types
- Fine-grained cytoarchitectural mapping of the DMN using a 3D reconstruction of human brain histology and non-linear manifold learning
- Comparison of DMN subregions based on the topography of their cytoarchtectural patterns (i.e. smooth or wavy)

**connectivity_and_cytoarchitecture.m**
- Relationship between connectivity and the principle axis of cytoarchitectural variation in the DMN, and how this varies depending on cortical type
- Evaluates the balance of extrinsic connectivity of the DMN, and other functional networks, with respect to cortical types
- This code generally works with any parcellated connectome. Data is made available to replicate our analyses, which focus on navigation efficiency derived from diffusion-weighted tractography as well as afferent and efferent connectivity derived from dynamic causal modelling of resting-state fMRI

### 🎆 Dependencies
The following open resources are necessary to run the code in this repository. Please note the code has been developed and tested using MATLAB2022a.
- BigBrainWarp (https://bigbrainwarp.readthedocs.io/en/latest/)
- micaopen (https://github.com/MICA-MNI/micaopen)
- Freesurfer, specifically the fsaverage5 subject directory (https://surfer.nmr.mgh.harvard.edu/)
- Gifti matlab (https://github.com/gllmflndn/gifti)
- BrainSpace (https://github.com/MICA-MNI/BrainSpace)
- Functional network atlas: Yeo2011_7Networks_N1000 (https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/)
- Von Economo atlas on fsaverage (http://www.dutchconnectomelab.nl/economo/)

 
### 📞 Contact
Feel free to get in touch if you have any questions: casey.paquola@gmail.com
