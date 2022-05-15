function e_dis = compute_external_force_distribution(j,q_new,q_old,qd_old,dt,time)
% the force will be exerted to the card at the first node in the first 1000
% points.
% At the 1001 point, the force will be canceled immediately.
global N l
Kd = 0.0001;
Kp = 0.005;
e_dis = zeros(4*N-1,1);
f = 3.7e-5;
if j< round(time/dt/100)
   e_dis(4*20-3:4*20-1) = [-f,0,0];
elseif j< round(time/dt/5)+1
   e_dis(1:3) = [0,0,4*f+Kp*(q_new(4*N-1)-q_new(3)-l+0.03)-Kd*qd_old(3)];

end


