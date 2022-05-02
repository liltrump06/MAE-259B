from parallel_transport import parallel_transport
from signedAngle import signedAngle
def computeReferenceTwist(u1, u2, t1, t2):
    ut = parallel_transport(u1, t1, t2)
    refTwist = signedAngle(ut, u2, t2)
    return refTwist