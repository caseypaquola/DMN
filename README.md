THE DEFAULT MODE NETWORK
----------------------------------------------

![Alt text](figures/figure_methods_overview.png?raw=true "Title")

If you've found your way to this GitHub, then I assume you've read our paper "The Unique Cytoarchiecture and Wiring of the Human Default Mode Network". Cool! If you would like dig deeper into the methods or reproduce aspects of the analyses, then please follow through these scripts.


### Key scripts
cytoarchitecture.m
- Evaluates the heterogeneity of the DMN according to cortical types
- Fine-grained cytoarchitectural mapping of the DMN using a 3D reconstruction of human brain histology and non-linear manifold learning
- Comparison of DMN subregions based on the topography of their cytoarchtectural patterns (i.e. smooth or wavy)

connectivity_and_cytoarchitecture.m
- Relationship between connectivity and the principle axis of cytoarchitectural variation in the DMN, and how this varies depending on cortical type
- Evaluates the balance of extrinsic connectivity of the DMN, and other functional networks, with respect to cortical types
- This code generally works with any parcellated connectome. Data is made available to replicate our analyses, which focus on navigation efficiency derived from diffusion-weighted tractography as well as afferent and efferent connectivity derived from dynamic causal modelling of resting-state fMRI

 
### Contact
Feel free to get in touch if you have any questions: casey.paquola@gmail.com
