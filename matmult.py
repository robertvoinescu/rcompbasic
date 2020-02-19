import numpy as np
import listconst as lst

def mult(aL,aO,aOV,bL,bO,bOV,N,guess,mto):
    Y=np.zeros([len(aL),len(bL)])
    # one electron integrals
    #########################
    # case: i=j
    # alpha strings
    for i in range(0,N):
        for I in range(0,len(aO.get(i))):
            a_key=aO[i][I]
            for b_key in bL:
                Y[a_key-1,b_key-1]+=guess[a_key-1,b_key-1]#*mto.get((i,i))
    # beta strings
    for i in range(0,N):
        for I in range(0,len(bO.get(i))):
            b_key=bO[i][I]
            for a_key in aL:
                Y[a_key-1,b_key-1]+=guess[a_key-1,b_key-1]#*mto.get((i,i))

    # case: i>j
    # alpha strings
    for i in range(0,N):
        for j in range(0,i):
            for I in range(0,len(aOV.get((i,j)))):
                a_key1 = aOV[(i,j)][I] 
                a_key2 = aOV[(j,i)][I] 
                s1 = -1 if a_key1<0 else 1
                s2 = -1 if a_key2<0 else 1
                a_key1=abs(a_key1)
                a_key2=abs(a_key2)
                for b_key in bL:
                    Y[a_key1-1,b_key-1]+=guess[a_key2-1,b_key-1]*s2
                    Y[a_key2-1,b_key-1]+=guess[a_key1-1,b_key-1]*s1
    
    for i in range(0,N):
        for j in range(0,i):
            for I in range(0,len(bOV.get((i,j)))):
                b_key1 = bOV[(i,j)][I] 
                b_key2 = bOV[(j,i)][I] 
                s1 = -1 if b_key1<0 else 1
                s2 = -1 if b_key2<0 else 1
                b_key1=abs(b_key1)
                b_key2=abs(b_key2)
                for b_key in bL:
                    Y[a_key-1,b_key1-1]+=guess[a_key-1,b_key2-1]*s2
                    Y[a_key-1,b_key2-1]+=guess[a_key-1,b_key1-1]*s1

    # two electron integrals
    #########################
    # diagonal terms
    for i in range(0,N):
        for I in range(0,len(aO.get(i))):
            a_key1=aO.get(i)[I]
            for J in range(0,len(bO.get(i))):
                bkey=bO.get(i)[J]
                Y[a_key-1,b_key-1]+=guess[a_key-1,bkey-1]
    
    # A^{alpha, beta} terms
    for j in range(0,N):
        for J in range(0,len(O.get(j))):
            key1=O.get(j)[J]
            for i in range(0,N):
                for k in range(0,i):
                    for I in range(0,len(OV.get((i,k)))):
                        key2=OV.get(i,k)[I]
                        key3=OV.get(k,i)[I]
                        s2=-1 if key2<0 else 1
                        s3=-1 if key3<0 else 1
                        Y[key1-1,key2]=

                        
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
    N=N_orbitals
    h={}
    for i in range(0,N):
        h[(i,i)]=1
    guess=np.ones([len(L),len(L)])
    print(mult(L,O,OV,L,O,OV,N,guess,h))
