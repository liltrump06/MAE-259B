function [F,J] = ComputeF_J_p1(N,q_new,q_old,qd_old,M,e_dis)
global dt dampratio dl E A I
g = [0,9.81]';
F = zeros(2*N,1);
J = zeros(2*N,2*N);
% stretching force and jacobian
for i = 1:N-1
    ind = 2*i-1:2*i+2;
    dF = gradEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2), dl, E*A);
    dJ = hessEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2), dl, E*A);
    F(ind) = F(ind)+dF;
    J(ind,ind) = J(ind,ind)+dJ;
end
% bending force and jacobian
for i = 1:N-2
    ind = 2*i-1:2*i+4;
    dF = gradEb(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),q_new(2*i+3),q_new(2*i+4),0,dl,E*I);
    dJ = hessEb(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),q_new(2*i+3),q_new(2*i+4),0,dl,E*I);
    F(ind) = F(ind)+dF;
    J(ind,ind) = J(ind,ind)+dJ;
end
% gravity and damping force
for m = 1:N 
    grav = M(2*m-1:2*m,2*m-1:2*m)*g;
    F(2*m-1:2*m) = F(2*m-1:2*m)+grav;
    damp = dampratio*(q_new(2*i-1:2*i)-q_old(2*i-1:2*i))/dt;
    F(2*m-1:2*m) = F(2*m-1:2*m)+damp;
end
% external force and inertia
F = F-e_dis;
F = F + M*(((q_new-q_old)/dt-qd_old)/dt);
% whole Jacobian
J = J + M*(1/dt)*(1/dt) + dampratio/dt * eye(size(J));