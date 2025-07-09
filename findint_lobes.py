"""
DESCRIPTION:
Given an NTD dimer and two chain ids, find and select the atoms on the first chain that are within 4.5 A of the second. 
PDBs are saved for whole interfaces and lobes individually.

USAGE:
findint pdbcode, chain1, chain2, reversechains, LLstart, LLend, pdbdir

PARAMS:
pdbcode (string of length 4)
		accession number in PDB of structure containing the NTD dimer
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
LLstart and LLend (numbers)
		Most of the lower lobe of iGluR NTDs and other PBP-related proteins is in one continuous segment with the remaining part being far from the dimer interface 
		so the upper lobe for the purposes of interface selection can be defined as all residues outside this segment. 
		If not given these default to 110 and 244, which Madhav used for 3hsy. 
pdbdir (string)
        directory with pdbs in it and where the results will be written
        defaults to current working directory

RETURNS:
6 selections and PDB files with the atoms in chain 1 within 4.5 A of chain 2 and vice versa as one interface and divided into upper and lower lobes by residue number.

EXAMPLE:
# find atoms in chain B of 3hsy within 4.5 A of chain A of 3hsy
findint 3hsy, B, A

AUTHOR:
James Krieger
"""

# change this for your system
import os
from pymol import cmd, stored

def findint(pdbcode, chain1='A', chain2='B', reversechains=-1,
            LLstart=110, LLend=244, pdbdir=os.getcwd()):
    # remove waters and hydrogens before starting
    cmd.remove("resn HOH")
    cmd.remove("hydro")

    if int(reversechains) == 0:
        dimer = pdbcode + "_" + chain1 + chain2
    elif int(reversechains) == 1:
        dimer = pdbcode + "_" + chain2 + chain1
    elif int(reversechains) == -1:
        dimer = pdbcode

    # set the name of the selection to return.
    rSelName1 = "IF_" + pdbcode + '_' + chain1
    ulSelName1 = 'UL_' + rSelName1
    llSelName1 = 'LL_' + rSelName1
    rSelName2 = "IF_" + pdbcode + '_' + chain2
    ulSelName2 = 'UL_' + rSelName2
    llSelName2 = 'LL_' + rSelName2

    # input checking
    if not checkParams(dimer):
        print("There was an error with a parameter.  Please see "
              "the above error message for how to fix it.")
        return None

    and_chain_str = ' and chain '
    open_br_str = '('
    close_br_str = ')'
    cmd.select(rSelName1, open_br_str + dimer + and_chain_str + chain1 + close_br_str +
               " within 4.5 of " + open_br_str + dimer + and_chain_str + chain2 + close_br_str)
    cmd.select(ulSelName1, rSelName1 + ' and not resi ' +
               str(LLstart) + '-' + str(LLend))
    cmd.select(llSelName1, rSelName1 + ' and resi ' +
               str(LLstart) + '-' + str(LLend))
    cmd.select(rSelName2, open_br_str + dimer + and_chain_str + chain2 + close_br_str +
               " within 4.5 of " + open_br_str + dimer + and_chain_str + chain1 + close_br_str)
    cmd.select(ulSelName2, rSelName2 + ' and not resi ' +
               str(LLstart) + '-' + str(LLend))
    cmd.select(llSelName2, rSelName2 + ' and resi ' +
               str(LLstart) + '-' + str(LLend))

    iter_str1 = "(byres %s) and n. ca"
    iter_str2 = "stored.reslist.append((resi,resn))"

    cmd.save(os.path.join(pdbdir, rSelName1 + ".pdb"), "{0}".format(rSelName1))
    f = open(os.path.join(pdbdir, rSelName1 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate(iter_str1 % rSelName1, iter_str2)
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, rSelName1 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(rSelName1, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + rSelName1 + "'"))

    cmd.save(os.path.join(pdbdir, ulSelName1 + ".pdb"), "{0}".format(ulSelName1))
    f = open(os.path.join(pdbdir, ulSelName1 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate(iter_str1 % ulSelName1, iter_str2)
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, ulSelName1 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(ulSelName1, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + ulSelName1 + "'"))

    cmd.save(os.path.join(pdbdir, llSelName1 + ".pdb"), "{0}".format(llSelName1))
    f = open(os.path.join(pdbdir, llSelName1 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate(iter_str1 % llSelName1, iter_str2)
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, llSelName1 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(llSelName1, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + llSelName1 + "'"))

    cmd.save(os.path.join(pdbdir, rSelName2 + ".pdb"), "{0}".format(rSelName2))
    f = open(os.path.join(pdbdir, rSelName2 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate(iter_str1 % rSelName2, iter_str2)
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, rSelName2 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(rSelName2, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + rSelName2 + "'"))

    cmd.save(os.path.join(pdbdir, ulSelName2 + ".pdb"), "{0}".format(ulSelName2))
    f = open(os.path.join(pdbdir, ulSelName2 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate(iter_str1 % ulSelName2, iter_str2)
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, ulSelName2 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(ulSelName2, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + ulSelName2 + "'"))

    cmd.save(os.path.join(pdbdir, llSelName2 + ".pdb"), "{0}".format(llSelName2))
    f = open(os.path.join(pdbdir, llSelName2 + "-residue-list.txt"), "w")
    stored.reslist = []
    cmd.iterate("(byres " + llSelName2 + ") and n. ca",
                "stored.reslist.append((resi,resn))")
    for i in range(len(stored.reslist)):
        f.write(str(stored.reslist[i]) + '\n')
    f.close()
    f = open(os.path.join(pdbdir, llSelName2 + "-atom-list.txt"), "w")
    stored.atomlist = []
    cmd.iterate(llSelName2, "stored.atomlist.append((resi,resn,name))")
    for i in range(len(stored.atomlist)):
        f.write(str(stored.atomlist[i]) + '\n')
    f.close()
    print(cmd.count_atoms("'" + llSelName2 + "'"))

    return rSelName1, ulSelName1, llSelName1, rSelName2, ulSelName2, llSelName2


cmd.extend("findint", findint)


def checkParams(dimer):
    """
    This is just a helper function for checking the user input
    """
    # check dimer
    if len(dimer) == 0 or not isinstance(dimer, str):
        print("Error: Please provide valid PyMOL object or selection name for dimer.")
        return False

    return True
