

## Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu)
# You should use this code at your own risk. Copy and redistribution is not
# permitted. Written permission is required.
import numpy as np
def computekappa(node0, node1, node2, m1e, m2e, m1f, m2f ):
#
# Inputs:
# node0: 1x3 vector - position of the node prior to the "turning" node
# node1: 1x3 vector - position of the "turning" node
# node2: 1x3 vector - position of the node after the "turning" node
#
# m1e: 1x3 vector - material director 1 of the edge prior to turning
# m2e: 1x3 vector - material director 2 of the edge prior to turning
# m1f: 1x3 vector - material director 1 of the edge after turning
# m2f: 1x3 vector - material director 2 of the edge after turning
#
# Outputs:
# kappa: 1x2 vector - curvature at the turning node
    t0 = (node1-node0) / np.linalg.norm(node1-node0)
    t1 = (node2-node1) / np.linalg.norm(node2-node1)
 
    kb = 2.0 * np.cross(t0.reshape(3), t1.reshape(3)) / (1.0 + np.inner(t0.T, t1.T))
    kappa = np.zeros([1, 2])
    kappa1 = 0.5 * np.inner( kb, m2e + m2f)
    kappa2 = -0.5 * np.inner( kb, m1e + m1f)
    kappa[:,0] = kappa1
    kappa[:,1] = kappa2
    return kappa
