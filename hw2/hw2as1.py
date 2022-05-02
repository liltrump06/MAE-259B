import numpy as np 
from computeReferenceTwist import computeReferenceTwist
## assignment1
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

twist1 = computeReferenceTwist(m11,m21,t1,t2)
twist2 = computeReferenceTwist(m21,m31,t2,t3)
print(twist1)
print(twist2)
