import numpy as np 
from crossMat import crossMat
def gradEs_hessEs(node0, node1,l_k, EA):
#
# Inputs:
# node0: 1x3 vector - position of the first node
# node1: 1x3 vector - position of the last node
#
# l_k: reference length (undeformed) of the edge
# EA: scalar - stretching stiffness - Young's modulus times area
#
# Outputs:
# dF: 6x1  vector - gradient of the stretching energy between node0 and node 1.
# dJ: 6x6 vector - hessian of the stretching energy between node0 and node 1.

## Gradient of Es
    edge = np.transpose(node1 - node0) #3x1 edge vector
    edgeLen = np.linalg.norm(edge)
    tangent = edge / edgeLen
    epsX = edgeLen/l_k - 1
    dF_unit = EA * tangent * epsX
    dF = np.zeros([6,1])
    dF[0:3] = - dF_unit
    dF[3:6] =   dF_unit
    ## Hessian of Es
    Id3 = np.eye(3)
    M = EA * ((1/l_k - 1/edgeLen) * Id3 +1/edgeLen * (edge*np.transpose(edge))/ np.float_power(edgeLen,2)) #Note edge * edge.T must be 3x3
    dJ = np.zeros([6,6])
    dJ[0:3, 0:3] = M
    dJ[3:6, 3:6] = M
    dJ[0:3, 3:6] = - M
    dJ[3:6, 0:3]= - M
    return dF,dJ


def gradEb_hessEb(node0, node1, node2,m1e, m2e, m1f, m2f,kappaBar, l_k, EI):
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
# kappaBar: 1x2 vector - natural curvature at the turning node
# l_k: voronoi length (undeformed) of the turning node
# EI: scalar - bending stiffness
#
# Outputs:
# dF: 11x1  vector - gradient of the bending energy at node1.
# dJ: 11x11 vector - hessian of the bending energy at node1.

## Computation of gradient of the two curvatures
    m1e =np.array([m1e])
    m2e =np.array([m2e])
    m1f =np.array([m1f])
    m2f =np.array([m2f])

    gradKappa = np.zeros([11,2])

    ee = node1 - node0
    ef = node2 - node1
 
    norm_e = np.linalg.norm(ee)
    norm_f = np.linalg.norm(ef)

    te = ee / norm_e
    tf = ef / norm_f

# Curvature binormal
    kb = 2.0 * np.cross(te.reshape(3), tf.reshape(3)) / (1.0 + np.dot(te, tf.T))

    chi = 1.0 + np.dot(te, tf.T)
    tilde_t = (te + tf) / chi
    tilde_d1 = (m1e + m1f) / chi
    tilde_d2 = (m2e + m2f) / chi

# Curvatures
    kappa1 = 0.5 * np.dot( kb, (m2e + m2f).T) # CHECKED
    kappa2 = -0.5 * np.dot( kb, (m1e + m1f).T) # CHECKED

# kappa1 = kappa(c, 1)
# kappa2 = kappa(c, 2)

    Dkappa1De = (1.0 / norm_e * (-kappa1 * tilde_t + np.cross(tf.reshape(3),tilde_d2.reshape(3))))
    Dkappa1Df = (1.0 / norm_f * (-kappa1 * tilde_t - np.cross(te.reshape(3),tilde_d2.reshape(3))))

    Dkappa2De = (1.0 / norm_e * (-kappa2 * tilde_t - np.cross(tf.reshape(3),tilde_d1.reshape(3))))
    Dkappa2Df = (1.0 / norm_f * (-kappa2 * tilde_t + np.cross(te.reshape(3),tilde_d1.reshape(3))))

    gradKappa[0:3, 0] = -Dkappa1De
    gradKappa[4:7, 0] = Dkappa1De - Dkappa1Df
    gradKappa[8:11, 0] = Dkappa1Df

    gradKappa[0:3, 1] = -Dkappa2De
    gradKappa[4:7, 1] = Dkappa2De - Dkappa2Df
    gradKappa[8:11, 1] = Dkappa2Df

    gradKappa[3, 0] = -0.5 * np.dot(kb, m1e.T)
    gradKappa[7, 0] = -0.5 * np.dot(kb, m1f.T)
    gradKappa[3, 1] = -0.5 * np.dot(kb, m2e.T)
    gradKappa[7, 1] = -0.5 * np.dot(kb, m2f.T)

    ## Computation of hessian of the two curvatures
    DDkappa1 = np.zeros([11, 11])
    DDkappa2 = np.zeros([11, 11])

    norm2_e = np.float_power(norm_e,2)
    norm2_f = np.float_power(norm_f,2)


    tt_o_tt = np.dot(tilde_t.T , tilde_t) # must be 3x3. tilde_t is 1x3
    tmp = np.array([np.cross(tf.reshape(3), tilde_d2.reshape(3))])
    tf_c_d2t_o_tt = np.dot(tmp.T , tilde_t) # must be 3x3
    tt_o_tf_c_d2t = tf_c_d2t_o_tt.T # must be 3x3
    kb_o_d2e = np.dot(kb.T , m2e) # must be 3x3
    d2e_o_kb = kb_o_d2e.T # must be 3x3

    Id3 = np.eye(3)
    D2kappa1De2 = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t)- kappa1/ (chi * norm2_e) * (Id3 - te.T*te) + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb)

    tmp = np.array([np.cross(te.reshape(3), tilde_d2.reshape(3))])
    te_c_d2t_o_tt = np.dot(tmp.T , tilde_t)
    tt_o_te_c_d2t = te_c_d2t_o_tt.T
    kb_o_d2f = kb.T * m2f
    d2f_o_kb = kb_o_d2f.T

    D2kappa1Df2 = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t)- kappa1 / (chi * norm2_f) * (Id3 - tf.T*tf)+ 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb)

    D2kappa1DeDf = -kappa1/(chi * norm_e * norm_f) * (Id3 + np.dot(te.T,tf))+ 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t - crossMat(tilde_d2))
    D2kappa1DfDe = D2kappa1DeDf.T

    tmp = np.array([np.cross(tf.reshape(3), tilde_d1.reshape(3))])
    tf_c_d1t_o_tt = np.dot(tmp.T,tilde_t) # must be 3x3
    tt_o_tf_c_d1t = tf_c_d1t_o_tt.T # must be 3x3
    kb_o_d1e = kb.T*m1e # must be 3x3
    d1e_o_kb = kb_o_d1e.T # must be 3x3

    D2kappa2De2= 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t)- kappa2 / (chi * norm2_e) * (Id3 - te.T*te)- 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb)

    tmp = np.array([np.cross(te.reshape(3), tilde_d1.reshape(3))])
    te_c_d1t_o_tt = np.dot(tmp.T,tilde_t) # must be 3x3
    tt_o_te_c_d1t = te_c_d1t_o_tt.T # must be 3x3
    kb_o_d1f = kb.T*m1f # must be 3x3
    d1f_o_kb =  kb_o_d1f.T # must be 3x3

    D2kappa2Df2 = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t)- kappa2 / (chi * norm2_f) * (Id3 - tf.T*tf)- 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb) # must be 3x3

    D2kappa2DeDf = -kappa2/(chi * norm_e * norm_f) * (Id3 + te.T*tf)+ 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1))# must be 3x3
    D2kappa2DfDe = D2kappa2DeDf.T # must be 3x3

    D2kappa1Dthetae2 = -0.5 * np.dot(kb, m2e.T)
    D2kappa1Dthetaf2 = -0.5 * np.dot(kb, m2f.T)
    D2kappa2Dthetae2 =  0.5 * np.dot(kb, m1e.T)
    D2kappa2Dthetaf2 =  0.5 * np.dot(kb, m1f.T)

    D2kappa1DeDthetae = 1.0 / norm_e * (0.5 * np.dot(kb, m1e.T) * tilde_t - 1.0 / chi * np.cross(tf.reshape(3), m1e.reshape(3)))
    D2kappa1DeDthetaf = 1.0 / norm_e * (0.5 * np.dot(kb, m1f.T) * tilde_t - 1.0 / chi * np.cross(tf.reshape(3), m1f.reshape(3)))
    D2kappa1DfDthetae = 1.0 / norm_f * (0.5 * np.dot(kb, m1e.T) * tilde_t + 1.0 / chi * np.cross(te.reshape(3), m1e.reshape(3)))
    D2kappa1DfDthetaf = 1.0 / norm_f * (0.5 * np.dot(kb, m1f.T) * tilde_t + 1.0 / chi * np.cross(te.reshape(3), m1f.reshape(3)))
    D2kappa2DeDthetae = 1.0 / norm_e * (0.5 * np.dot(kb, m2e.T) * tilde_t - 1.0 / chi * np.cross(tf.reshape(3), m2e.reshape(3)))
    D2kappa2DeDthetaf = 1.0 / norm_e * (0.5 * np.dot(kb, m2f.T) * tilde_t - 1.0 / chi * np.cross(tf.reshape(3), m2f.reshape(3)))
    D2kappa2DfDthetae = 1.0 / norm_f * (0.5 * np.dot(kb, m2e.T) * tilde_t + 1.0 / chi * np.cross(te.reshape(3), m2e.reshape(3)))
    D2kappa2DfDthetaf = 1.0 / norm_f * (0.5 * np.dot(kb, m2f.T) * tilde_t + 1.0 / chi * np.cross(te.reshape(3), m2f.reshape(3)))

# Curvature terms
    DDkappa1[0:3, 0:3]  =   D2kappa1De2
    DDkappa1[0:3, 4:7]  = - D2kappa1De2 + D2kappa1DeDf
    DDkappa1[0:3, 8:11] =               - D2kappa1DeDf
    DDkappa1[4:7, 0:3]  = - D2kappa1De2                + D2kappa1DfDe
    DDkappa1[4:7, 4:7]  =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2
    DDkappa1[4:7, 8:11]=                 D2kappa1DeDf                - D2kappa1Df2
    DDkappa1[8:11, 0:3]  =                              - D2kappa1DfDe
    DDkappa1[8:11, 4:7]  =                                D2kappa1DfDe - D2kappa1Df2
    DDkappa1[8:11, 8:11] =                                               D2kappa1Df2

# Twist terms
    DDkappa1[3, 3]     =   D2kappa1Dthetae2
    DDkappa1[7, 7]     =   D2kappa1Dthetaf2

# Curvature-twist coupled terms
    DDkappa1[0:3, 3]   = - D2kappa1DeDthetae
    DDkappa1[4:7, 3]   =   D2kappa1DeDthetae - D2kappa1DfDthetae
    DDkappa1[8:11,3]   =                       D2kappa1DfDthetae
    DDkappa1[3, 0:3]   =   np.transpose(DDkappa1[0:3, 3])
    DDkappa1[3, 4:7]   =   np.transpose(DDkappa1[4:7, 3])
    DDkappa1[3, 8:11]  =   np.transpose(DDkappa1[8:11,3])

# Curvature-twist coupled terms
    DDkappa1[0:3, 7]   = - D2kappa1DeDthetaf
    DDkappa1[4:7, 7]   =   D2kappa1DeDthetaf - D2kappa1DfDthetaf
    DDkappa1[8:11, 7]  =                       D2kappa1DfDthetaf
    DDkappa1[7, 0:3]   =   np.transpose(DDkappa1[0:3, 7])
    DDkappa1[7, 4:7]   =   np.transpose(DDkappa1[4:7, 7])
    DDkappa1[7, 8:11]  =   np.transpose(DDkappa1[8:11,7])

# Curvature terms
    DDkappa2[0:3, 0:3] =   D2kappa2De2
    DDkappa2[0:3, 4:7] = - D2kappa2De2 + D2kappa2DeDf
    DDkappa2[0:3, 8:11] =               - D2kappa2DeDf
    DDkappa2[4:7, 0:3] = - D2kappa2De2                + D2kappa2DfDe
    DDkappa2[4:7, 4:7] =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2
    DDkappa2[4:7, 8:11]=                 D2kappa2DeDf                - D2kappa2Df2
    DDkappa2[8:11, 0:3]=                              - D2kappa2DfDe
    DDkappa2[8:11, 4:7]=                                D2kappa2DfDe - D2kappa2Df2
    DDkappa2[8:11, 8:11]=                                               D2kappa2Df2

# Twist terms
    DDkappa2[3, 3]     = D2kappa2Dthetae2
    DDkappa2[7, 7]     = D2kappa2Dthetaf2

# Curvature-twist coupled terms
    DDkappa2[0:3, 3]   = - D2kappa2DeDthetae
    DDkappa2[4:7, 3]   =   D2kappa2DeDthetae - D2kappa2DfDthetae
    DDkappa2[8:11,3]   =                       D2kappa2DfDthetae
    DDkappa2[3, 0:3]   =   np.transpose(DDkappa2[0:3, 3])
    DDkappa2[3, 4:7]   =   np.transpose(DDkappa2[4:7, 3])
    DDkappa2[3, 8:11]  =   np.transpose(DDkappa2[8:11,3])

# Curvature-twist coupled terms
    DDkappa2[0:3, 7]   = - D2kappa2DeDthetaf
    DDkappa2[4:7, 7]  =   D2kappa2DeDthetaf - D2kappa2DfDthetaf
    DDkappa2[8:11,7]   =                       D2kappa2DfDthetaf
    DDkappa2[7, 0:3]   =   np.transpose(DDkappa2[0:3, 7])
    DDkappa2[7, 4:7]   =   np.transpose(DDkappa2[4:7, 7])
    DDkappa2[7,8:11]   =   np.transpose(DDkappa2[8:11,7])
    
## Gradient of Eb
    EIMat = np.eye(2)*EI
    kappaVector = np.array([[float(kappa1),float(kappa2)]])
    dkappaVector = kappaVector - kappaBar
    dF = np.dot(np.dot(gradKappa , EIMat) , dkappaVector.T) / l_k

## Hessian of Eb
    dJ = 1.0 / l_k * np.dot(np.dot(gradKappa , EIMat) , np.transpose(gradKappa))
    temp = np.array(1.0 / l_k * np.dot(dkappaVector , EIMat))
    dJ = dJ + (temp[0,0] * DDkappa1 + temp[0,1] * DDkappa2)
    return dF,dJ

def gradEt_hessEt(node0, node1, node2,theta_e, theta_f, refTwist,l_k, GJ):
#
# Inputs:
# node0: 1x3 vector - position of the node prior to the "twisting" node
# node1: 1x3 vector - position of the "twisting" node
# node2: 1x3 vector - position of the node after the "twisting" node
#
# theta_e: scalar - twist angle of the first edge
# theta_f: scalar - twist angle of the second (last) edge
#
# l_k: voronoi length (undeformed) of the turning node
# GJ: scalar - twisting stiffness
#
# Outputs:
# dF: 11x1  vector - gradient of the twisting energy at node1.
# dJ: 11x11 vector - hessian of the twisting energy at node1.
## Computation of gradient of the twist
    gradTwist = np.zeros([11,1])
    ee = node1 - node0
    ef = node2 - node1
    norm_e = np.linalg.norm(ee)
    norm_f = np.linalg.norm(ef)
    norm2_e = np.float_power(norm_e,2)
    norm2_f = np.float_power(norm_f,2)
    te = ee / norm_e
    tf = ef / norm_f
# Curvature binormal
    kb = 2.0 * np.cross(te, tf) / (1.0 + np.dot(te, tf.T))
    gradTwist[0:3] = -0.5 / norm_e * kb.T
    gradTwist[8:11] = 0.5 / norm_f * kb.T
    gradTwist[4:7] = -(gradTwist[0:3]+gradTwist[8:11])
    gradTwist[3] = -1
    gradTwist[7] = 1
## Computation of hessian of twist
    DDtwist = np.zeros([11, 11])
    chi = 1.0 + np.dot(te, tf.T)
    tilde_t = (te + tf) / chi
    D2mDe2 = -0.25 / norm2_e * ( kb.T * (te + tilde_t)+ (te + tilde_t).T * kb)
    D2mDf2 = -0.25 / norm2_f  * ( kb.T * (tf + tilde_t) + (tf + tilde_t).T * kb )
    D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - kb.T * tilde_t )
    D2mDfDe = D2mDeDf.T
    DDtwist[0:3,0:3] = D2mDe2
    DDtwist[0:3, 4:7] = -D2mDe2 + D2mDeDf
    DDtwist[4:7, 0:3] = -D2mDe2 + D2mDfDe
    DDtwist[4:7, 4:7] = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2
    DDtwist[0:3, 8:11] = -D2mDeDf
    DDtwist[8:11, 0:3] = -D2mDfDe
    DDtwist[8:11, 4:7] = D2mDfDe - D2mDf2
    DDtwist[4:7,8:11] = D2mDeDf - D2mDf2
    DDtwist[8:11,8:11] = D2mDf2
    
## Gradient of Et
    integratedTwist = theta_f - theta_e + refTwist
    dF = GJ/l_k * integratedTwist * gradTwist
## Hessian of Eb
    dJ = GJ/l_k * (integratedTwist * DDtwist + gradTwist*gradTwist.T)
    return dF,dJ 