import numpy as np
import listconst as lst

n_electrons=2
N_orbitals=4

L=lst.Lgen(N_orbitals,n_electrons)
O=lst.Ogen(L,N_orbitals)
OV=lst.OVgen(L,N_orbitals)
OO=lst.OOgen(L,N_orbitals)
OOV=lst.OOVgen(L,N_orbitals)
OOVV=lst.OOVVgen(L,N_orbitals)

