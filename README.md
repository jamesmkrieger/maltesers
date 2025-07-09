# maltesers
A set of python pymol scripts for carrying out the dimer interface analysis originally done in  Sukumaran et al., EMBO J 2011

These scripts now allow this analysis to be done in PyMOL and easily generate the final images.

1. Run PyMOL and load a structure with all these scripts in a directory it can see, e.g. run PyMOL from inside the maltesers directory

    pymol example/4gpa_AB.pdb

2. Add the maltesers command to PyMOL by selecting the menu option File -> Run Script... and selecting maltesers.py

3. Run the maltesers command on a structure loaded into PyMOL as follows. The pdbdir should be relative to the maltesers directory or an absolute path

    maltesers 4gpa_AB, pdbdir=example

After running the maltesers command in PyMOL, you end up with two images, which can then be combined or overlaid, e.g. in Powerpoint
