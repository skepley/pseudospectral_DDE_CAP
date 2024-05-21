MATLAB code associated with our paper on validation of unstable manifolds for pseudospectral approximations of delay differential equations which can be found [here](
https://arxiv.org/abs/2405.07727). Below are instructions for reproducing the Wright's equation example from that paper. 

1. Install MATLAB. This code was produced and tested using MATLAB 2019a. 
2. Clone the IMP library available [here](https://github.com/skepley/IMP) to a folder in your MATLAB path.
3. Be sure to have INTLAB installed and available in your MATLAB path. This code was produced and tested using INTLAB v.9.
4. Clone the pseudospectral_DDE_CAP repo (this is the page you are at). Make sure it is available in your MATLAB path.
5. Type "startintlab" in MATLAB and press Enter.
6. Open the file named "wright_eqn_ode_validation.m". Change parameters at the top of this script (if desired) and then run it. The code takes ~6-10 hours and its recommended to save the output once it completes. Alternatively, you can load "ode_data_final.m" to obtain the output of this file for the default parameters chosen in the paper.
7. Open the file named "wright_eqn_dde_validation.m" and run this file to perform the validation of the DDE using the approximate manifold of the PSA.
