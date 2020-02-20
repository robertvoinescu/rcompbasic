import numpy as np
import listconst as lst

def mult(L,O,OV,OO,OOV,N,X,mto):
    Y=np.zeros([len(L),len(L)])
    # one electron integrals
    #########################
    # case: i=j
    for i in range(0,N):
        for I in range(0,0 if O.get(i) == None else len(O.get(i))):
            o_key=O[i][I]
            for b_key in L:
                Y[o_key-1,b_key-1]+=X[o_key-1,b_key-1]#*mto.get((i,i))
            for a_key in L:
                Y[a_key-1,o_key-1]+=X[a_key-1,o_key-1]#*mto.get((i,i))

    # case: i>j
    # alpha strings
    for i in range(0,N):
        for j in range(0,i):
            for I in range(0,0 if OV.get((i,j)) == None else len(OV.get((i,j)))):
                ov_key1 = OV[(i,j)][I] 
                ov_key2 = OV[(j,i)][I] 
                s1 = -1 if ov_key1<0 else 1
                s2 = -1 if ov_key2<0 else 1
                ov_key1=abs(ov_key1)
                ov_key2=abs(ov_key2)
                for b_key in L:
                    Y[ov_key1-1,b_key-1]+=X[ov_key2-1,b_key-1]*s2
                    Y[ov_key2-1,b_key-1]+=X[ov_key1-1,b_key-1]*s1
                for a_key in L:
                    Y[a_key-1,ov_key1-1]+=X[a_key-1,ov_key2-1]*s2
                    Y[a_key-1,ov_key2-1]+=X[a_key-1,ov_key1-1]*s1

    # two electron integrals
    #########################
    # diagonal terms
    for i in range(0,N):
        for I in range(0,0 if O.get(i) == None else len(O.get(i))):
            o_key1=O.get(i)[I]
            for J in range(0,0 if O.get(i) == None else len(O.get(i))):
                o_key2=O.get(i)[J]
                Y[o_key1-1,o_key2-1]+=X[o_key1-1,o_key2-1]
    for i in range(0,N):
        for j in range(0,i):
            for I in range(0,0 if O.get(i) == None else len(O.get(i))):
                o_key1=O.get(i)[I]
                for J in range(0,0 if O.get(j) == None else len(O.get(j))):
                    o_key2=O.get(j)[J]
                    Y[o_key1-1,o_key2-1]+=X[o_key1-1,o_key2-1]
                    Y[o_key2-1,o_key1-1]+=X[o_key2-1,o_key1-1]

    
    # A^{alpha, beta} terms
    for j in range(0,N):
        for J in range(0,0 if O.get(j) == None else len(O.get(j))):
            o_key=O.get(j)[J]
            for i in range(0,N):
                for k in range(0,i):
                    for I in range(0,0 if OV.get((i,k)) == None else len(OV.get((i,k)))):
                        ov_key1=OV[(i,k)][I]
                        ov_key2=OV[(k,i)][I]
                        s1=-1 if ov_key1<0 else 1
                        s2=-1 if ov_key2<0 else 1
                        ov_key1=abs(ov_key1)
                        ov_key2=abs(ov_key2)
                        Y[o_key-1,ov_key1-1]+=X[o_key-1,ov_key2-1]*s2
                        Y[o_key-1,ov_key2-1]+=X[o_key-1,ov_key1-1]*s1
                        Y[ov_key1-1,o_key-1]+=X[ov_key2-1,o_key-1]*s2
                        Y[ov_key2-1,o_key-1]+=X[ov_key1-1,o_key-1]*s1

    #B^{alpha, beta}
    for i in range(0,N):
        for k in range(0,i):
            for l in range(0,N):
                for j in range(0,l):
                    for I in range(0,0 if OV.get((i,k)) == None else len(OV.get((i,k)))):
                        for J in range(0,0 if OV.get((j,l)) == None else len(OV.get((j,l)))):
                            ov_key1=OV[(i,k)][I]
                            ov_key2=OV[(k,i)][I]
                            ov_key3=OV[(j,l)][J]
                            ov_key4=OV[(l,j)][J]
                            s1=-1 if ov_key1<0 else 1
                            s2=-1 if ov_key3<0 else 1
                            ov_key1=abs(ov_key1)
                            ov_key2=abs(ov_key2)
                            ov_key3=abs(ov_key3)
                            ov_key4=abs(ov_key4)
                            Y[ov_key2-1,ov_key4-1]+=X[ov_key1-1,ov_key3-1]*s1*s2
                            Y[ov_key2-1,ov_key3-1]+=X[ov_key1-1,ov_key4-1]*s1*s2
                            Y[ov_key1-1,ov_key4-1]+=X[ov_key2-1,ov_key3-1]*s1*s2
                            Y[ov_key1-1,ov_key3-1]+=X[ov_key2-1,ov_key4-1]*s1*s2

    # Insert some off diagnol terms here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # alpha-alpha and beta-beta cases                  
    # diagonal terms which involve the lists OO(i,j) (i>j)
    for i in range(0,N):
        for j in range(0,j):
            V=1
            for I in range(0,0 if OO.get((i,j)) == None else len(OO.get((i,j)))):
                oo_key=OO[(i,j)][I]
                for b_key in L:
                    Y[oo_key-1,b_key-1]+=X[oo_key-1,b_key-1]*V
                for a_key in L:
                    Y[a_key-1,oo_key-1]+=+X[a_key-1,oo_key-1]*V
                    
    # three index terms with OOV(i,j,k)
    for j in range(0,N):
        for k in range(0,N):
            for i in [x for x in range(0,N) if x != j and x != k]:
                V=1
                for I in range(0,0 if OOV.get((i,j)) == None else len(OOV.get((i,j)))):
                    oov_key1=OOV[(i,j,k)][I]
                    oov_key2=OOV[(i,k,j)][I]
                    s1=-1 if oov_key1<0 else 1
                    s2=-1 if oov_key3<0 else 1
                    oov_key1=abs(oov_key1)
                    oov_key2=abs(oov_key2)
                    for b_key in L:
                        Y[oov_key1-1,b_key-1]+=X[oov_key2-1,b_key-1]*s2*V
                        Y[oov_key2-1,b_key-1]+=X[oov_key1-1,b_key-1]*s1*V
                    for a_key in L:
                        Y[a_key-1,oov_key1-1]+=X[a_key-1,oov_key2-1]*s2*V
                        Y[a_key-1,oov_key2-1]+=X[a_key-1,oov_key1-1]*s1*V

    # four index terms associated with OOVV(i,j,k,l) i>j, k>l
    for i in range(0,N):
        for j in range(0,i):
            for k in range(0,i):
                for l in range(0,min(k,j)):
                    if i!=j and i!=k and i!=l and j!=k and j!=l and k!=l:
                        V=1
                        for I in range(0,0 if OOVV.get((i,j)) == None else len(OOVV.get((i,j)))):
                            oovv_key1=OOVV[(i,j,k,l)][I]
                            oovv_key2=OOVV[(k,l,i,j)][I]
                            s1=-1 if oovv_key1<0 else 1
                            s2=-1 if oovv_key3<0 else 1
                            oov_key1=abs(oovv_key1)
                            oov_key2=abs(oovv_key2)
                            for b_key in L:
                                y[oovv_key1-1,b_key-1]+=X[oovv_key2-1,b_key-1]*s2*V
                                y[oovv_key2-1,b_key-1]+=X[oovv_key1-1,b_key-1]*s1*V
                            for a_key in L:
                                y[a_key-1,oovv_key1-1]+=X[a_key-1,oovv_key2-1]*s2*V
                                y[a_key-1,oovv_key2-1]+=X[a_key-1,oovv_key1-1]*s1*V

    return Y

                





if __name__=='__main__':
    import time
    n_electrons=2
    N_orbitals=8
    
    start_time0=time.time()
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
    start_time1=time.time()
    guess=np.ones([len(L),len(L)])
    print('Our vector is of length {}'.format(len(L)**2))
    print(mult(L,O,OV,OO,OOV,N,guess,h))
    print('The lists where created in {} and the matrix multiplication took {}.'.format(time.time()-start_time0,time.time()-start_time1))
