
from signedAngle import signedAngle
from parallel_transport import parallel_transport
from computekappa import computekappa
from computeF_J import computeF_J
import math
import numpy as np

N = 50
q = np.zeros([4*N-1,1],np.float64)
qd = np.zeros([4*N-1,1],np.float64)
kappabarlist = np.zeros([N-2,2],np.float64)
dt = 0.01
l = 0.2
Rn = 0.02
r0 = 0.001
dtheta = l/Rn*1/(N-1)
rou = 1000
E = 1e7
G = E/3
g = np.array([[0,0,-9.81]],np.float64)
J = np.pi*np.float_power(r0,4)/2
I = np.pi*np.float_power(r0,4)/4
A = np.pi*np.float_power(r0,2)
EA = E*A
EI = E*I
GJ = G*J
M = np.eye(4*N-1)
for i in range(0,N):
    if i == N-1:
        xk = np.array([[Rn* math.cos(i*dtheta), Rn*math.sin(i*dtheta),0]])
        q[4*i:4*i+3] = xk.T
        M[4*i:4*i+3,4*i:4*i+3] = np.eye(3)*A*l*rou/(N)
    else:
        xk = np.array([[Rn* math.cos(i*dtheta), Rn*math.sin(i*dtheta),0]])
        q[4*i:4*i+3] = xk.T
        q[4*i+3] = 0
        M[4*i:4*i+3,4*i:4*i+3] = np.eye(3)*A*l*rou/(N)
        M[4*i+3,4*i+3] = A*l*rou/(N-1)*r0*r0/2
l_k = (np.linalg.norm(q[8:11]-q[4:7])+np.linalg.norm(q[4:7]-q[0:3]))/2

##parameter prepare
u1 = np.zeros([N-1,3],np.float64)
u2 = np.zeros([N-1,3],np.float64)
tk_old = np.zeros([N-1,3],np.float64)
tk_new = np.zeros([N-1,3],np.float64)
a1_old = np.zeros([N-1,3],np.float64)
a1_new = np.zeros([N-1,3],np.float64)
a2_old = np.zeros([N-1,3],np.float64)
a2_new = np.zeros([N-1,3],np.float64)
q_new = q
q_old = q
m1 = np.zeros([N-1,3],np.float64)
m2 = np.zeros([N-1,3],np.float64)
mkref = np.zeros([N-2,1],np.float64)
err = 100
free_index = slice(7,4*N)
print(free_index)
t1 = (q[4:7] - q[0:3])/np.linalg.norm((q[4:7] - q[0:3]))
u11 = np.array([[-t1[1,0],t1[0,0],0]])
a11 = u11
time = 5
q_whole = np.zeros([4*N-1,int(time/dt)],np.float64)
qd_whole = np.zeros([4*N-1,int(time/dt)],np.float64)
for i in range(0,N-2):
    kappabarlist[i,:] = [0.2047929168960332,0]

for j in range(0,500):
    print('j=',j)
    err = 100
    while err >1e-4:
        for m in range(0,N-1):
            tk1 =(q_new[4*m+4:4*m+7] - q_new[4*m:4*m+3])/np.linalg.norm((q_new[4*m+4:4*m+7] - q_new[4*m:4*m+3]))
            tk1 = tk1.reshape(3)
            if j == 0:
                if m == 0:
                    u1[m,:] = u11
                else:
                    u1[m,:] = parallel_transport(u1[m-1,:],tk_new[m-1,:],tk1)
                a1_new[m,:] = u1[m,:]
            else:
                if m == 0:
                    u1[m,:] = u11
                else:
                    u1[m,:] = parallel_transport(u1[m-1,:],tk_new[m-1,:],tk1)
                a1_new[m,:] = parallel_transport(a1_old[m,:],tk_old[m,:],tk1)
            a2_new[m,:] = np.cross(tk1,a1_new[m,:])
            u2[m,:] = np.cross(tk1,u1[m,:])
            tk_new[m,:] = tk1
            m1[m,:] = math.cos(q_new[4*m+3])*a1_new[m,:]+math.sin(q_new[4*m+3])*a2_new[m,:]
            m2[m,:] = math.cos(q_new[4*m+3])*a2_new[m,:]-math.sin(q_new[4*m+3])*a1_new[m,:]

        for k in range(0,N-2):
            mkref[k,0] = signedAngle(parallel_transport(a1_new[k,:],tk_new[k,:],tk_new[k+1,:]),a1_new[k+1,:],tk_new[k,:])
        [F,Jaco] = computeF_J(N,q_new,l_k,EA,m1,m2,GJ,EI,kappabarlist,mkref,g,M)
        F = F+np.dot(M,((q_new-q_old)/dt-qd)/dt)
        Jaco = Jaco + M*(1/dt)*(1/dt) 
        F_free = F[free_index]
        J_free = Jaco[free_index,free_index]
        Jinv = np.linalg.inv(J_free)
        deltaq = np.matmul(Jinv,F_free)

        q_new[free_index] = q_new[free_index] - deltaq 
        err = np.sum(np.abs(F[free_index]))
        print('err=',err)  
    q_whole[:,j] = q_new.reshape(4*N-1)
    qd = (q_new-q_old)*(1/dt)
    q_old = q_new
    a1_old = a1_new
    a2_old = a2_new
    tk_old = tk_new

print(a1_old)




