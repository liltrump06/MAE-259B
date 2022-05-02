import numpy as np
from Elastic_Rod import gradEt_hessEt, gradEs_hessEs,gradEb_hessEb
node0 = np.array([[0.1,0.1,0.1]])
node1 = np.array([[0.2,0.2,0.1]])
node2 = np.array([[0.3,0.2,0.1]])
l_k =0.03
EA = 1e8
m1e = np.array([0,0,1])
m2e = np.array([0,1,0])
m1f = np.array([0,0.707,-0.707])
m2f = np.array([0,-0.707,0.707])
kappaBar = np.array([0.2,0.1])
EI = 1e8
theta_e = 0.2
theta_f = 0.4
GJ = 1e8
reftwist = 0.3
[dF,dJ] = gradEs_hessEs(node0, node1,l_k, EA)
print(dF)
print(dJ)
[dF,dJ] = gradEb_hessEb(node0, node1, node2,m1e, m2e, m1f, m2f,kappaBar, l_k, EI)
print(dF)
print(dJ)
[dF,dJ] = gradEt_hessEt(node0, node1, node2,theta_e, theta_f, reftwist,l_k, GJ)
print(dF)
print(dJ)