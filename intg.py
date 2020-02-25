from functools import reduce
import numpy
from pyscf import ao2mo
from pyscf import __config__
def h2id(i,j,k,l):
    if i >= j:
        ij = i * (i-1) // 2 + j-1
    else:
        ij = j * (j-1) // 2 + i-1
    if k >= l:
        kl = k * (k-1) // 2 + l-1
    else:
        kl = l * (l-1) // 2 + k-1
    if ij >= kl:
        index=ij*(ij+1)//2+kl 
    else:
        index=kl*(kl+1)//2+ij
    return index
if __name__=='__main__':
    from pyscf import gto, scf, fci, tools
    mol = gto.mole(atom = 'h 0 0 0; f 0 0 1.1', basis= '6-31g', symmetry=true)
    myhf = scf.rhf(mol)
    myhf.kernel()
    cisolver = fci.fci(mol, myhf.mo_coeff)
    tools.fcidump.from_scf(myhf,'temp.txt')
    print(tools.fcidump.read('temp.txt')['h2'][0])        
    print(tools.fcidump.read('temp.txt')['h2'][h2id(0,8,0,2)])        
