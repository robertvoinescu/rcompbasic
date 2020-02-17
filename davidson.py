#!/bin/python
import numpy as np
import time


n = 1200            					
tol = 1e-8				
mmax = n//2


# INPUT
sparsity = 0.0001
A = np.zeros((n,n))
for i in range(0,n):
    A[i,i] = i + 1 
A = A + sparsity*np.random.randn(n,n) 
A = (A.T + A)/2 


k = 8			 
eig = 4			
t = np.eye(n,k)		
V = np.zeros((n,n))	
I = np.eye(n)		


start_davidson = time.time()

for m in range(k,mmax,k):
    if m <= k:
        for j in range(0,k):
            V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
        theta_old = 1 
    elif m > k:
        theta_old = theta[:eig]
    V,R = np.linalg.qr(V)
    T = np.dot(V[:,:(m+1)].T,np.dot(A,V[:,:(m+1)]))
    THETA,S = np.linalg.eig(T)
    idx = THETA.argsort()
    theta = THETA[idx]
    s = S[:,idx]
    for j in range(0,k):
        w = np.dot((A - theta[j]*I),np.dot(V[:,:(m+1)],s[:,j])) 
        q = w/(theta[j]-A[j,j])
        V[:,(m+j+1)] = q
    norm = np.linalg.norm(theta[:eig] - theta_old)
    if norm < tol:
        break

end_davidson = time.time()


print("I found the following eigenv", theta[:eig]," in ",
    end_davidson - start_davidson, "seconds")

