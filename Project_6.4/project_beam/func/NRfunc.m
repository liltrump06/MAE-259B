function NRfunc(q_new,velop)
q_p = q_new;
while err >1e-5
    e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time);
    [F,J] = ComputeF_J(N,q_new,M,e_dis,dl,E,A,I);
    F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
    Jini = M*(1/dt)*(1/dt);
    J = J + Jini;
    deltaX = J(free_index,free_index) \ F(free_index);
    q_new(free_index) = q_new(free_index) - deltaX;
    err = sum(abs(F(free_index)));
end