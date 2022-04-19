function J = ComputeJ(q_new,l,E,A,I,M,C,dti)
JEs1 = hessEs(q_new(1),q_new(2),q_new(3),q_new(4),l/2,E*A);
JEs1_full = [JEs1,zeros(4,2);zeros(2,6)];
JEs2 = hessEs(q_new(3),q_new(4),q_new(5),q_new(6),l/2,E*A);
JEs2_full = [zeros(2,6);zeros(4,2),JEs2];
JEb = hessEb(q_new(1),q_new(2),q_new(3),q_new(4),q_new(5),q_new(6),0,l/2,E*I);
JE = JEs1_full+JEb+JEs2_full;
for k = 1:6
    for j = 1:6
        if k == j
            J_ine(k,j) = M(k,k)/dti^2;
            J_vis(k,j) = C(k,k)/dti;
        else 
            J_ine(k,j) = 0;
            J_vis(k,j) = 0;
        end
    end
end
J =JE+ J_ine + J_vis;