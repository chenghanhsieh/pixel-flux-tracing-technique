# Pixel Flux-tracing Technique (PTF)
Last Updated: April 20 2023 by Cheng-Han Hsieh (cheng-han.hsieh@yale.edu) 
All rights reserved by Cheng-Han Hsieh @Yale.

This is the main code for the Pixel Flux-tracing Technique developed by Cheng-Han Hsieh et al. 2023. in prep.

#Prerequisite
Astropy , Numpy, python 2.7

# What will this page contains? 
This folder contains functions to compute 2D outflow mass, momentum, energy map using the Pixel Flux-tracing Technique developed by Hsieh et al. 2023. in prep.


In this folder, I include an example for HH212. The HH212 CO J = 3-2 data is observed by the Atacama Large Millimeter/submillimeter Array (ALMA).

To run the script, please download the HH212 fits cube: 

[https://drive.google.com/file/d/1zwwOQr1Lw03PPmb4pfTQCfMvqPwM8WGt/view?usp=sharing
](https://drive.google.com/file/d/1ZPmo8jLalSR5k7365lUoLolU3BRbIuSE/view?usp=sharing)

1. First run the density script to get a H2 column density cube. 

python CO_density.py

The code will generate Den_HH212_CO32.fits. 

2. Then run the Pixel Flux-tracing Technique code. This code will generate the 2D instanteous mass, momentum, and energy rate maps for HH212. 

python dEdMdP_maps.py


If you use this method for your research, please cite: Hsieh et al. 2023. 
