function e_dis = compute_external_force(j,q_new,q_old,qd_old,dt,time)
% the force will be exerted to the card at the first node in the first 1000
% points.
% At the 1001 point, the force will be canceled immediately.
global N l Rc
Kd = 0.01;
Kp = 0.05;
e_dis = zeros(2*N,1);
f = 9e-3;
if j< round(time/dt/100)
   e_dis(2*10-1:2*10) = [-f,0];
else
   e_dis(1:2) = [0,4*f+Kp*(q_new(2*N)-q_new(2)-l+0.03)-Kd*qd_old(2)];
end
