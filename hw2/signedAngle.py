import numpy as np
def signedAngle(u,v,n):
# "angle" is signed angle from vector "u" to vector "v" with axis "n"
    w = np.cross(u,v)
    angle = np.arctan2(np.linalg.norm(w),np.dot(u,v))
    if (np.dot(n,w)<0):
        angle = -angle
    return angle