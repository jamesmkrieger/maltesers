''' module for use in pymol to calculate upper lobe and lower lobe contact numbers
 syntax in pymol is as follows
 CBF pdbcode,chain1,chain2,lobe,reversechains

reversechains can take 0 or 1 or -1. 0 is the default
1 means use the chains in the reverse order for the object name e.g. CBF 3saj,C,A,UL,1 would be used if you had an object called 3saj_AC
0 means don't reverse the order e.g. CBF 3saj,A,C,UL is fine
-1 means don't use them at all - maybe your pdb structure only has one dimer so you haven't split it into objects with chain labels e.g. CBF 3hsy,A,B,UL,-1
'''
from pymol import *

def CBF(pdbcode,chain1,chain2,lobe,reversechains=0):
	fi = open('Documents/pdbs/' + lobe + '_IF_' + pdbcode + '_' + chain1 + '.pdb','r')
	for line in fi:
		if line.find('ATOM') == 0:
			resid = line.split()[5]
			atomName = line.split()[2]
			if int(reversechains) == 0:
				fo = open('Documents/pdbs/' + lobe + '_IF_' + pdbcode + '_' + chain1 + chain2 + '_' + chain1 + '_CBF','a')
				fo.write(line[:61] + str(cmd.count_atoms(pdbcode + "_" + chain1 + chain2 + " and chain " + chain2 + " within 4.5 of /" + pdbcode + "_" + chain1 + chain2 + "//" + chain1 + "/" + resid + "/" + atomName)).rjust(5) + line[66:])
				fo.close()
			elif int(reversechains) == 1:
                                fo = open('Documents/pdbs/' + lobe + '_IF_' + pdbcode + '_' + chain2 + chain1 + '_' + chain1 + '_CBF','a')
                                fo.write(line[:61] + str(cmd.count_atoms(pdbcode + "_" + chain2 + chain1 + " and chain " + chain2 + " within 4.5 of /" + pdbcode + "_" + chain2 + chain1 + "//" + chain1 + "/" + resid + "/" + atomName)).rjust(5) + line[66:])
				fo.close()
			elif int(reversechains) == -1:
				fo = open('Documents/pdbs/' + lobe + '_IF_' + pdbcode + '_' + chain1 + '_CBF','a')
                                fo.write(line[:61] + str(cmd.count_atoms(pdbcode + " and chain " + chain2 + " within 4.5 of /" + pdbcode + "//" + chain1 + "/" + resid + "/" + atomName)).rjust(5) + line[66:])
				fo.close()
	fi.close()

cmd.extend("CBF",CBF)