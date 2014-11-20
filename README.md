The scattering library
======================

Routines for computing the analytical solutions to acoustic scattering problems involving prolate spheroids, oblate spheroids, and disks.

The software in this library is described in the following paper: http://scitation.aip.org/content/asa/journal/jasa/136/6/10.1121/1.4901318.

## Using the scattering library ##

1. Clone this repository to your local machine.
Alternatively, you can download this repository as a ZIP by clicking on the link to the right.
This will give you the entire repository, including not only the scattering library, but also the spheroidal library.

2. Open up MATLAB, and run `sandbox.m`.
This will set up all of the paths.

3. Run `generate_all_test_figures.m`.
This repository comes with the precomputed spheroidal wave functions necessary to generate these test figures.
It will generate four figures, saving each one as a PDF in the `images` directory.

4. To run `generate_figures.m`, you'll need to download more precomputed spheroidal wave functions at https://dl.dropboxusercontent.com/u/83368249/saved.zip.
This ZIP is approximately 800 MB, so make sure you have enough disk space.
Unzip the contents of saved.zip into the `saved` directory.
Once done, run `generate_all_figures.m`, which will generate 11 figures, saving each one as a PDF in the `images` directory.

5. For those interested in lower-frequency stuff, you can download https://dl.dropboxusercontent.com/u/83368249/saved_small_c.zip, which contains precomputed spheroidal wave functions for c = 0.1, 0.5, 1.0, and 5.0.

6. You can use the the spheroidal library to compute the spheroidal wave functions for values of c, m, and/or n that aren't included in either of these ZIPs.

## License ##

The scattering library is Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami, and is released under the BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause).

The spheroidal library
======================

The spheroidal library includes routines for computing the spheroidal wave functions, and is located in the `spheroidal` directory.

Third-party software
====================

This repository contains four libraries from the MATLAB Central File Exchange.
They are:

1. `export_fig` by Oliver Woodford (http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)

2. `cubehelix` by Stephen Cobeldick (http://www.mathworks.com/matlabcentral/fileexchange/43700-cubehelix-colormaps)

3. `legendflex` by Kelly Kearney (http://www.mathworks.com/matlabcentral/fileexchange/31092-legendflex--a-more-flexible-legend)

4. `spaceplots` by Aditya (http://www.mathworks.com/matlabcentral/fileexchange/35464-spaceplots)

These four libraries were released under their own licenses, which can be found along with their code in their respective directories in this repository.
