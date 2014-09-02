scattering
==========

Routines for computing the analytical solutions to acoustic scattering problems involving prolate spheroids, oblate spheroids, and disks.

## Using the scattering library ##

1. Clone this repository to your local machine.
Alternatively, you can download this repository as a ZIP by clicking on the link to the right.

2. Open up MATLAB, and run `sandbox.m`.
This will set up all of the paths.

3. Run `generate_all_test_figures.m`.
This repository comes with the precomputed spheroidal wave functions necessary to generate these test figures.
It will generate four figures, saving each one as a PDF in the `images` directory.

4. To run `generate_figures.m`, you'll need to download more precomputed spheroidal wave functions at https://dl.dropboxusercontent.com/u/83368249/saved.zip.
This ZIP is approximately 800 MB, so make sure you have enough disk space.
Unzip the contents of saved.zip into the `saved` directory.
Once done, run `generate_figures.m`, which will generate 11 figures, saving each one as a PDF in the `images` directory.
