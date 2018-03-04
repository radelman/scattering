The spheroidal library
======================

Routines for computing the spheroidal wave functions.

The software in this library is described in the following paper: http://arxiv.org/abs/1408.0074.

## Using the spheroidal library ##

1. Clone this repository to your local machine.
Alternatively, you can download this repository as a ZIP by navigating to the root directory of this repository and clicking on the link to the right.
This will give you the entire repository, including not only the spheroidal library, but also the scattering library.

2. To use the spheroidal library, you will need the binaries for `pro_sphwv` and `obl_sphwv`.
The source code for these is in the `sphwv` directory.
In order to compile them on your own, however, you must download the GNU MPFR library (http://www.mpfr.org).
Since this can be tedious to do, we have compiled them for you, and provide the binaries at the following two links.
For the Linux binaries (specifically, RHEL 5.5), download https://www.dropbox.com/s/y19qu4ifjcy6bjh/sphwv_linux.zip?dl=0.
For the Windows binaries, download https://www.dropbox.com/s/a6dswosf5e4fji9/sphwv_windows.zip?dl=0.
Whichever platform you're on, the corresponding ZIP will contain two executables.
Copy both of them into the `sphwv` directory.

2. Open up MATLAB, enter the `spheroidal` directory, and run `sandbox.m`.
This will set up all of the paths.

3. Run `check_sphwv.m`.
This will check that `pro_sphwv` and `obl_sphwv` are working, and should only take a few seconds.

4. Run `generate_data_for_all_test_figures.m`, followed by `generate_all_test_figures.m`.
This will compute some test data, and then generate two test figures, saving each one as a PDF in the `images` directory.
While very accurate, the use of arbitrary precision arithmetic causes the software to be slow.
On one of our machines, generating the data takes around an hour.

5. Run `generate_data_for_all_figures.m`, followed by `generate_all_figures.m`.
This will compute some data, and then generate 14 figures, saving each one as a PDF in the `images` directory.
Again, generating this data is slow, and will take around two hours.

6. The scripts, `generate_tasks.m`, `run_tasks.m`, and `run_tasks.pl`, allow you to easily generate and run scripts for computing the spheroidal wave functions for different values of c, m, and n.  For more information, see the comments in these scripts.

## License ##

The spheroidal library is Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami, and is released under the BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause).
