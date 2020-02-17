import numpy as np
import listconst as lst

def mult(aL,aO,bL,bO,alpha,beta,guess):
    Y={}
    for i in range(0,len(alpha)):
        for I in range(0,len(aO.get(i))):
            a_key=aO[i][I]
            for b_key in bL:
                Y[(a_key,b_key)]=guess.get(a_key,b_key)*h.get((i,i))
    return Y

                





if __name__=='__main__':
    n_electrons=2
    N_orbitals=4

    L=lst.Lgen(N_orbitals,n_electrons)
    O=lst.Ogen(L,N_orbitals)
    OV=lst.OVgen(L,N_orbitals)
    OO=lst.OOgen(L,N_orbitals)
    OOV=lst.OOVgen(L,N_orbitals)
    OOVV=lst.OOVVgen(L,N_orbitals)
    alpha=[1]*N_orbitals
    beta=[1]*N_orbitals
    h={}
    from pyscf import gto, scf, fci
    mol = gto.Mole(atom = 'H 0 0 0; F 0 0 1.1', basis= '6-31g', symmetry=True)
    myhf = scf.RHF(mol)
    myhf.kernel()
    cisolver = fci.FCI(mol, myhf.mo_coeff)
    cisolver.kernel()
    for i in range(0,len(alpha)):
        h[(i,i)]=1
    print(mult(L,O,L,O,alpha,beta,h))
