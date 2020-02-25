import numpy as np

def proj(u,v):
    proj=u*np.dot(v,u)/np.dot(u,u)
    return proj

def gram_schmidt(M,g):
    for gs in range(g,M.shape[0]):
        for bs in range(0,gs):
            M[:,gs]-=proj(M[:,bs],M[:,gs])
        M[:,gs]=M[:,gs]/np.linalg.norm(M[:,gs])
    return M

if __name__=='__main__':
    g=2
    M=np.array([[1,0,2,1],[0,1,15,0],[0,0,1,0],[0,0,0,1]],dtype=float)            
    print('The following list of vectors:\n',M)
    print('is decomposed along the first {} columns.'.format(g))
    print(gram_schmidt(M,g))        

        
    
