function q_new = Compute_q_new_exp(q_old,qd_old,M,C,W,l,E,A,I,dti)

dEs1 = gradEs(q_old(1),q_old(2),q_old(3),q_old(4),l/2,E*A);
dEs2 = gradEs(q_old(3),q_old(4),q_old(5),q_old(6),l/2,E*A);
dEb = gradEb(q_old(1),q_old(2),q_old(3),q_old(4),q_old(5),q_old(6),0,l/2,E*I);
dE = [dEs1;0;0]+[0;0;dEs2]+dEb;
q_new = q_old + dti*(inv(M)*dti*(-dE-C*qd_old+W)+qd_old);
end