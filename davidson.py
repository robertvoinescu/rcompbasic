#!/bin/pythonj
import listconst as lst
import numpy as np
import time
import math
import matmult as mt
from pyscf import tools

def nCr(n,r):
    if n<r:
        comb=0
    else:
        f = math.factorial
        comb=f(n) / f(r) / f(n-r)
    return int(comb)

def dot(H,V):
    X=np.zeros([H.shape[0],V.shape[1]]) 
    for i in range(0,V.shape[1]):
        X[:,i]=np.dot(H,V[:,i])
    return X

def dotH(V,L,O,OV,OO,OOV,OOVV,fcidump):
    n=nCr(fcidump['NORB'],fcidump['NELEC']//2)
    X=np.zeros([n**2,V.shape[1]]) 
    for i in range(0,V.shape[1]):
        X[:,i]=np.reshape(mt.mult(L,O,OV,OO,OOV,OOVV,np.reshape(V[:,i],(n,n)),fcidump),(n**2))
    return X
tol = 1e-8				


filename='temp.txt'
fcidump=tools.fcidump.read(filename)
NELEC=int(fcidump['NELEC'])
NORB=int(fcidump['NORB'])
n=nCr(NORB,NELEC//2)**2

filename='temp.txt'
fcidump=tools.fcidump.read(filename)
NELEC=fcidump['NELEC']
NORB=fcidump['NORB']
print(NELEC,NORB)

start_time0=time.time()
L=lst.Lgen(NORB,NELEC//2)
O=lst.Ogen(L,NORB)
OV=lst.OVgen(L,NORB)
OO=lst.OOgen(L,NORB)
OOV=lst.OOVgen(L,NORB)
OOVV=lst.OOVVgen(L,NORB)
print(n)

# INPUT
sparsity = 0.0001
A = np.zeros((n,n))
for i in range(0,n):
    A[i,i] = i + 1 
A = A + sparsity*np.random.randn(n,n) 
A = (A.T + A)/2 


mmx=8
k = 8					# number of initial guess vectors 
eig = 4					# number of eignvalues to solve 
t = np.eye(n,k)			# set of k unit vectors as guess
V = np.zeros((n,k*mmx))		# array of zeros to hold guess vec

dotH(V,L,O,OV,OO,OOV,OOVV,fcidump)
# Begin block Davidson routine

start_davidson = time.time()

for m in range(k,(mmx-1)*k):
    if m <= k:
        for j in range(0,k):
            V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
        theta_old = 1 
    elif m > k:
        theta_old = theta[:eig]
    V,R = np.linalg.qr(V)
    sigma=dotH(V[:,:(m+1)],L,O,OV,OO,OOV,OOVV,fcidump)
    T = np.dot(V[:,:(m+1)].T,sigma)
    THETA,S = np.linalg.eig(T)
    idx = THETA.argsort()
    theta = THETA[idx]
    s = S[:,idx]
    for j in range(0,k):
        w = np.dot(sigma,s[:,j])-theta[j]*np.dot(V[:,:(m+1)],s[:,j])
        q = w/(theta[j]-sigma[j,j])
        V[:,(m+j+1)] = q
    norm = np.linalg.norm(theta[:eig] - theta_old)
    if norm < tol:
        break
    if m == (mmx-1)*k-1 or m==n-10:
        print('ERROR: storage vector not large enough. results are not within tol = ',tol)
        print('TOL = ',norm)
        break
end_davidson = time.time()

print("Eigenvalues: ", theta[:eig]," found in ",
    end_davidson - start_davidson, "seconds")

