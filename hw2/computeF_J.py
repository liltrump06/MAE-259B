import numpy as np
from Elastic_Rod import gradEb_hessEb, gradEs_hessEs, gradEt_hessEt

def computeF_J(N,q,l_k,EA,m1,m2,GJ,EI,kappalist,mkref,g,M):
    F = np.zeros([4*N-1,1],np.float64)
    Jaco = np.zeros([4*N-1,4*N-1],np.float64)
    for k in range(0,N-1):
        ind1 = slice(4*k,4*k+3)
        ind2 = slice(4*k+4,4*k+7)
        [dF,dJ] = gradEs_hessEs(q[4*k:4*k+3].T,q[4*k+4:4*k+7].T,l_k,EA)
        F[ind1] = F[ind1]+dF[0:3]
        F[ind2] = F[ind2]+dF[3:6]
        Jaco[ind1,ind1] = Jaco[ind1,ind1]+dJ[0:3,0:3]
        Jaco[ind1,ind2] = Jaco[ind1,ind2]+dJ[0:3,3:6]      
        Jaco[ind2,ind1] = Jaco[ind2,ind1]+dJ[3:6,0:3]
        Jaco[ind2,ind2] = Jaco[ind2,ind2]+dJ[3:6,3:6]
    for k in range(0,N-2):
        ind = slice(4*k,4*k+11)
        [dF,dJ] = gradEb_hessEb(q[4*k:4*k+3].T,q[4*k+4:4*k+7].T,\
            q[4*k+8:4*k+11].T,m1[k,:],m2[k,:],m1[k+1,:],m2[k+1,:],kappalist[k,:],l_k,EI)
        F[ind] = F[ind]+dF
        Jaco[ind,ind] = Jaco[ind,ind]+dJ
        [dF,dJ] = gradEt_hessEt(q[4*k:4*k+3].T,q[4*k+4:4*k+7].T,\
            q[4*k+8:4*k+11].T,q[4*k+3],q[4*k+7],mkref[k,:],l_k,GJ)
        F[ind] = F[ind]+dF
        Jaco[ind,ind] = Jaco[ind,ind]+dJ
    for k in range(0,N):
        F[4*k:4*k+3] = F[4*k:4*k+3]- np.matmul(M[4*k:4*k+3,4*k:4*k+3],g.T)
    return F,Jaco

