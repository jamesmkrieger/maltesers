''' module for use in pymol to calculate upper lobe and lower lobe contact numbers
 syntax in pymol is as follows
 LD pdbcode, chain1, chain2, lobe, reverse, pdbdir

reversechains can take 0 or 1 or -1. 0 is the default
1 means use the chains in the reverse order for the object name e.g. LD 3saj,C,A,UL,1 would be used if you had an object called 3saj_AC
0 means don't reverse the order e.g. LD 3saj,A,C,UL is fine if you have an object called 3saj_AC
-1 means don't use them at all - maybe your pdb structure only has one dimer so you haven't split it into objects with chain labels e.g. LD 3hsy,A,B,UL,-1
'''
from pymol import cmd
import os


def LD(pdbcode, chain1, chain2, lobe, reversechains=0,
       pdbdir=os.getcwd()):
    within_12_str = " within 12 of /"
    totalCount = 0.
    i = 0.
    fi = open(os.path.join(pdbdir, lobe + '_IF_' +
              pdbcode + '_' + chain1 + '.pdb'), 'r')
    fo = open(os.path.join(pdbdir, lobe + '_IF_' +
              pdbcode + '_' + chain1 + '_LD.txt'), 'w')
    for line in fi:
        if line.find('ATOM') == 0:
            resid = line.split()[5]
            atomName = line.split()[2]
            if int(reversechains) == 0:
                currCount = cmd.count_atoms(lobe + "_IF_" + pdbcode + "_" + chain1 + within_12_str +
                                            pdbcode + "_" + chain1 + chain2 + "//" + chain1 + "/" + resid + "/" + atomName)
            elif int(reversechains) == 1:
                currCount = cmd.count_atoms(lobe + "_IF_" + pdbcode + "_" + chain1 + within_12_str +
                                            pdbcode + "_" + chain2 + chain1 + "//" + chain1 + "/" + resid + "/" + atomName)
            elif int(reversechains) == -1:
                currCount = cmd.count_atoms(lobe + "_IF_" + pdbcode + "_" + chain1 +
                                            within_12_str + pdbcode + "//" + chain1 + "/" + resid + "/" + atomName)
            fo.write(line[:61] + str(currCount).rjust(5) + line[66:])
            totalCount += float(currCount)
            i += 1.
    fi.close()

    LDresult = totalCount / i
    print(LDresult)
    fo.write(str(LDresult))
    fo.close()


cmd.extend("LD", LD)
