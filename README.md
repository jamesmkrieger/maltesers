# maltesers
A set of python pymol scripts for carrying out the dimer interface analysis originally done in  Sukumaran et al., EMBO J 2011

These scripts now allow this analysis to be done in PyMOL and easily generate the final images.

1. Run PyMOL and load a structure with all these scripts in a directory it can see, e.g. run PyMOL from inside the maltesers directory

    pymol example/4gpa_AB.pdb

2. Add the maltesers command to PyMOL by selecting the menu option File -> Run Script... and selecting maltesers.py

3. Run the maltesers command on a structure loaded into PyMOL as follows. The pdbdir should be relative to the maltesers directory or an absolute path

    maltesers 4gpa_AB, pdbdir=example

After running the maltesers command in PyMOL, you end up with two images, which can then be combined with the following python commands or in an image editor:

    # combine the images for chain 1 to make the maltesers image
    spheres_1 = cv2.imread(pdbcode + midString1 + '_' + chain1 + '_spheres.png',1)
    surf_1 = cv2.imread(pdbcode + midString1 + '_' + chain1 + '_white_surface.png',1)
    malt_1 = cv2.addWeighted(spheres_1,0.5,surf_1,0.5,0)
    cv2.imwrite(pdbcode + midString1 + '_' + chain1 + '_maltesers.png',malt_1)

    # combine the images for chain 2 to make the maltesers image
    spheres_2 = cv2.imread(pdbcode + midString2 + '_' + chain2 + '_spheres.png',1)
    surf_2 = cv2.imread(pdbcode + midString2 + '_' + chain2 + '_white_surface.png',1)
    malt_2 = cv2.addWeighted(spheres_2,0.5,surf_2,0.5,0)
    cv2.imwrite(pdbcode + midString2 + '_' + chain2 + '_maltesers.png',malt_2)
