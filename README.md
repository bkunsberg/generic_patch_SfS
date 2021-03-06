# generic_patch_SfS
Matlab Code to Generate Generic Height Patches from a Polynomial Shaded Patch


Code implementing the algorithm described in Part II of 
Tensors, Differential Geometry and Statistical Shading Analysis
Daniel N. Holtmann-Rice, Benjamin S. Kunsberg, Steven W. Zucker
https://link.springer.com/article/10.1007/s10851-018-0815-z


Code written by Benjamin Kunsberg and tested under Matlab 2016, 2017, 2018.

To run, first unzip the Tensor Toolbox contained in tensor_toolbox.zip.  (Found here: https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html)


Then, run an example via the script: 'run_examples.m'

You can find more examples by commenting out the different lines in 'buildTestImg.m' or loading in your own image via Matlab's imread.  Note that the algorithm will only attempt to model the image in a central patch defined by the Taylor order you choose. 


