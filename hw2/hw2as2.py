import numpy as np 
from computeReferenceTwist import computeReferenceTwist
from parallel_transport import parallel_transport
from signedAngle import signedAngle
## assignment2
x1 = np.array([0,0,0]).T
x2 = np.array([0.5,0,0]).T
x3 = np.array([0.75,0.25,0]).T
x4 = np.array([0.75,0.5,0.25]).T
m11 = np.array([0,0,1]).T
m21 = np.array([0,0,1]).T
m31 = np.array([0,-np.sqrt(2)/2,np.sqrt(2)/2]).T
t1 = (x2-x1)/np.linalg.norm(x2-x1)
t2 = (x3-x2)/np.linalg.norm(x3-x2)
t3 = (x4-x3)/np.linalg.norm(x4-x3)
u11 = np.array([0,0,1]).T
u21 = parallel_transport(u11,t1,t2)
u31 = parallel_transport(u21,t2,t3)
ang1 = signedAngle(u11,m11,t1)
ang2 = signedAngle(u21,m21,t2)
ang3 = signedAngle(u31,m31,t3)
tau2 = ang2-ang1
tau3 = ang3-ang2
print('u31=',u31)
print('tau2=',tau2)
print('tau3=',tau3)
a = range(0,3)+1
print(a)
