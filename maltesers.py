"""
DESCRIPTION:
Given an NTD dimer and two chain IDs, do the maltesers analysis of Sukumaran et al., EMBO J 2011:
1. find the interface atoms (within 4.5 A) in both lobes using the findint function from findint_lobes.py
2. calculate the numbers of contacts and make an image showing them coloured red to blue by 
colouring by B-factor using the CBF function from CBF.py and also make an image of the interface surface
3. print out the local contact density as defined by Bahadur et al., J Mol Biol 2004 
using the LD function from LD.py

The process would need to be completed by combining the images together with opencv module cv2.

USAGE:
maltesers pdbcode, chain1, chain2, reversechains, pdbdir

PARAMS:
pdbcode (string)
                name of the object in pymol containing the NTD dimer
chain1 (string of length 1)
                chain id 1.
                defaults to A
chain2 (string of length 1)
                chain id 2
                defaults to B
reversechains (0, 1 or -1)
                reversechains can take 0 or 1 or -1. -1 is the default
                1 means use the chains in the reverse order for the object name e.g. 3saj,C,A,UL,1 would be used if you had an object called 3saj_AC
                0 means don't reverse the order e.g. 3saj,A,C,UL is fine if you have an object called 3saj_AC
                -1 means don't use them at all - maybe your pdb structure only has one dimer so you haven't split it into objects with chain labels e.g. 3hsy,A,B,UL,-1
pdbdir (string)
                directory with pdbs in it and where the results will be written
                defaults to current working directory

RETURNS:
A number of selections and PDB files as well as 4 png files

EXAMPLE:
# Do the maltesers analysis on dimer CD from 3o21
maltesers 3o21, A, B

AUTHOR:
James Krieger
"""
from LD import LD
from CBF import CBF
from findint_lobes import findint
import os
from pymol import cmd

def maltesers(pdbcode, chain1='A', chain2='B', reversechains=-1,
              pdbdir=os.getcwd()):
    
    if not os.path.exists(pdbdir):
        os.mkdir(pdbdir)

    cmd.disable()
    cmd.set('auto_zoom', 0)
    cmd.set('ray_shadows', 'off')

    if int(reversechains) == 0:
        midString1 = '_' + chain1 + chain2
        midString2 = '_' + chain2 + chain1
    elif int(reversechains) == 1:
        midString1 = '_' + chain2 + chain1
        midString2 = '_' + chain1 + chain2
    else:
        midString1 = midString2 = ''

    # use 3hsy as an orientation reference
    view_str = '''\
    -0.344349116,   -0.499494404,    0.794939816,\
     0.360954791,    0.711206675,    0.603238463,\
    -0.866680741,    0.494661897,   -0.064608492,\
     0.000000000,    0.000000000, -250.839691162,\
    15.650138855,   -2.656803131,  -12.363032341,\
   197.763916016,  303.915466309,  -20.000000000 '''
    if '3hsy' not in cmd.get_object_list():
        cmd.fetch('3hsy')
    cmd.set_view(view_str)
    cmd.align(pdbcode + midString1, '3hsy')
    cmd.disable('3hsy')

    # create the white surface image for chain 1
    cmd.bg_color('white')
    cmd.disable(pdbcode + midString1)
    cmd.create(pdbcode + midString1 + '_' + chain1,
               pdbcode + midString1 + ' and chain ' + chain1)
    cmd.show_as('surface', pdbcode + midString1 + '_' + chain1)
    cmd.color('white', pdbcode + midString1 + '_' + chain1)
    cmd.ray(1200, 1200)
    cmd.png(os.path.join(pdbdir, pdbcode + midString1 + '_' + chain1 + '_white_surface'))

    # create the white surface image for chain 2
    cmd.disable(pdbcode + midString1 + '_' + chain1)
    cmd.create(pdbcode + midString2 + '_' + chain2,
               pdbcode + midString2 + ' and chain ' + chain2)
    cmd.show_as('surface', pdbcode + midString2 + '_' + chain2)
    cmd.color('white', pdbcode + midString2 + '_' + chain2)
    cmd.set_view(view_str)
    cmd.turn('y', 180)
    cmd.ray(1200, 1200)
    cmd.png(os.path.join(pdbdir, pdbcode + midString2 + '_' + chain2 + '_white_surface'))

    findint(pdbcode, chain1, chain2, reversechains, pdbdir=pdbdir)

    cbd_pdb_str = '_CBF.pdb'

    # make the coloured spheres image for chain 1
    cmd.disable(pdbcode + midString2 + '_' + chain2)
    CBF(pdbcode, chain1, chain2, 'UL', reversechains, pdbdir=pdbdir)
    CBF(pdbcode, chain1, chain2, 'LL', reversechains, pdbdir=pdbdir)
    cmd.load(os.path.join(pdbdir, 'UL_IF_' + pdbcode +
             midString1 + '_' + chain1 + cbd_pdb_str))
    cmd.show_as('spheres', 'UL_IF_' + pdbcode +
                midString1 + '_' + chain1 + '_CBF')
    cmd.spectrum('b', 'blue_red', 'UL_IF_' + pdbcode +
                 midString1 + '_' + chain1 + '_CBF', 1, 7)
    if os.path.isfile(os.path.join(pdbdir, 'LL_IF_' + pdbcode + 
                                   midString1 + '_' + chain1 + cbd_pdb_str)):
        cmd.load(os.path.join(pdbdir, 'LL_IF_' + pdbcode +
                 midString1 + '_' + chain1 + cbd_pdb_str))
        cmd.show_as('spheres', 'LL_IF_' + pdbcode +
                    midString1 + '_' + chain1 + '_CBF')
        cmd.spectrum('b', 'blue_red', 'LL_IF_' + pdbcode +
                     midString1 + '_' + chain1 + '_CBF', 1, 7)
    cmd.set_view(view_str)
    cmd.ray(1200, 1200)
    cmd.png(os.path.join(pdbdir, pdbcode + midString1 + '_' + chain1 + '_spheres'))

    # make the coloured spheres image for chain 2
    cmd.disable('UL_IF_' + pdbcode + midString1 + '_' + chain1 + '_CBF')
    CBF(pdbcode, chain2, chain1, 'UL', reversechains, pdbdir=pdbdir)
    cmd.load(os.path.join(pdbdir, 'UL_IF_' + pdbcode +
             midString2 + '_' + chain2 + cbd_pdb_str))
    cmd.show_as('spheres', 'UL_IF_' + pdbcode +
                midString2 + '_' + chain2 + '_CBF')
    cmd.spectrum('b', 'blue_red', 'UL_IF_' + pdbcode +
                 midString2 + '_' + chain2 + '_CBF', 1, 7)
    if os.path.isfile(os.path.join(pdbdir, 'LL_IF_' + pdbcode + 
                                   midString1 + '_' + chain1 + cbd_pdb_str)):
        cmd.disable('LL_IF_' + pdbcode + midString1 + '_' + chain1 + '_CBF')
        CBF(pdbcode, chain2, chain1, 'LL', reversechains, pdbdir=pdbdir)
        cmd.load(os.path.join(pdbdir, 'LL_IF_' + pdbcode +
                 midString2 + '_' + chain2 + cbd_pdb_str))
        cmd.show_as('spheres', 'LL_IF_' + pdbcode +
                    midString2 + '_' + chain2 + '_CBF')
        cmd.spectrum('b', 'blue_red', 'LL_IF_' + pdbcode +
                     midString2 + '_' + chain2 + '_CBF', 1, 7)
    cmd.set_view(view_str)
    cmd.turn('y', 180)
    cmd.ray(1200, 1200)
    cmd.png(os.path.join(pdbdir, pdbcode + midString2 + '_' + chain2 + '_spheres'))

    # show the LD values
    ld_txt_str = '_LD.txt'
    LD(pdbcode, chain1, chain2, 'UL', reversechains, pdbdir=pdbdir)
    fi = open(os.path.join(pdbdir, 'UL_IF_' + pdbcode +
              midString1 + '_' + chain1 + ld_txt_str), 'r')
    lines = fi.readlines()
    fi.close()
    UL_LD_1 = float(lines[-1])

    LD(pdbcode, chain2, chain1, 'UL', reversechains, pdbdir=pdbdir)
    fi = open(os.path.join(pdbdir, 'UL_IF_' + pdbcode +
              midString2 + '_' + chain2 + ld_txt_str), 'r')
    lines = fi.readlines()
    fi.close()
    UL_LD_2 = float(lines[-1])

    print((UL_LD_1 + UL_LD_2)/2.)

    LD(pdbcode, chain1, chain2, 'LL', reversechains, pdbdir=pdbdir)
    fi = open(os.path.join(pdbdir, 'LL_IF_' + pdbcode +
              midString1 + '_' + chain1 + ld_txt_str), 'r')
    lines = fi.readlines()
    fi.close()
    LL_LD_1 = float(lines[-1])

    LD(pdbcode, chain2, chain1, 'LL', reversechains, pdbdir=pdbdir)
    fi = open(os.path.join(pdbdir, 'LL_IF_' + pdbcode +
              midString2 + '_' + chain2 + ld_txt_str), 'r')
    lines = fi.readlines()
    fi.close()
    LL_LD_2 = float(lines[-1])

    print((LL_LD_1 + LL_LD_2)/2.)

    cmd.delete('3hsy')


cmd.extend('maltesers', maltesers)
