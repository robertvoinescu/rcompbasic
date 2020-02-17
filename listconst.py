import math
import itertools


def nCr(n,r):
    if n<r:
        comb=0
    else:
        f = math.factorial
        comb=f(n) / f(r) / f(n-r)
    return comb 

def Add(sigma):
    lex=1
    occk=0
    for k in range(0,len(sigma)):
        if sigma[k]==1:
            occk+=1
        lex+=sigma[k]*nCr(k,occk)
    return int(lex)

def Lgen(N,n):
    L={}
    for bits in itertools.combinations(range(N),n):
        sigmab=[0]*N
        for bit in bits:
            sigmab[bit]=1
        L[Add(sigmab)]=sigmab
    return L

def Ogen(L,N):
    O={}
    for i in range(0,N):
        register=[]
        for key in L:
            if L[key][i]==1:
                register+=[key]
        if register:
            O[(i)]=register
    return O

def OVgen(L,N):
    def sign(i,j):
        if (i-j)%2==0:
            pm=1
        else:
            pm=-1
        return pm
    OV={}
    for i in range(0, N):
        for j in range(0, N): 
            pm=sign(i,j)
            register=[]
            for key in L:
                if [L[key][i]==1, L[key][j]]==[1,0]:
                    register+=[pm*key]
            if register:
                OV[(i,j)]=register
    return OV

def OOgen(L,N):
    OO={}
    for i in range(0,N):
        for j in range(i+1,N): 
            register=[]
            for key in L:
                if [L[key][i]==1, L[key][j]]==[1,1]:
                    register+=[key]
            if register:
                OO[(i,j)]=register
    return OO

def OOVgen(L,N):
    def sign(j,k,L,key):
        begin=min(j,k)+1
        end=max(j,k)
        if sum(L[key][begin:end])%2==0:
            pm=1
        else:
            pm=-1
        return pm

    OOV={}
    for i in range(0,N):
        for j in range(i+1, N): 
            for k in range(0,N):
                register=[]
                for key in L:
                    if [L[key][i]==1, L[key][j], L[key][k]]==[1,1,0]:
                        register+=[sign(j,k,L,key)*key]
                if register:
                    OOV[(i,j,k)]=register
    return OOV

def OOVVgen(L,N):
    def sign(i,j,k,l,L,key):
        b1=min(k,i)
        b2=min(l,j)
        e1=max(k,i)
        e2=max(l,j)
        if sum(L[key][b1:e1])%2==0 and  sum(L[key][b2:e2])%2==0:
            pm=1
        else:
            pm=-1
        return pm

    OOVV={}
    for i in range(0,N):
        for j in range(i+1, N): 
            for k in range(0, N):
                for l in range(k+1, N):
                    register=[]
                    for key in L:
                        if [L[key][i]==1, L[key][j], L[key][k], L[key][l]]==[1,1,0,0]:
                            register+=[sign(i,j,k,l,L,key)*key]
                    if register:
                        OOVV[(i,j,k,l)]=register
    return OOVV

if __name__ == '__main__':
    N=4
    n=2
    L=Lgen(N,n)
    print('Test code for {} orbitals with {} electrons:'.format(N,n))
    print('Our main list is: ', L)
    print('The O list is: {}'.format(Ogen(L,N)))
    print('The OV list is: {}'.format(OVgen(L,N)))
    print('The OO list is: {}'.format(OOgen(L,N)))
    print('The OOV list is: {}'.format(OOVgen(L,N)))
    print('The OOVV list is: {}'.format(OOVVgen(L,N)))

