from pyscf import gto, scf, fci, tools
mol = gto.Mole(atom = 'H 0 0 0; F 0 0 1.1', basis= '6-31g', symmetry=True)
myhf = scf.RHF(mol)
myhf.kernel()
cisolver = fci.FCI(mol, myhf.mo_coeff)
tools.fcidump.from_scf(myhf,'temp.txt')
