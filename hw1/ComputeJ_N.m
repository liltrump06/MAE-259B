function J = ComputeJ_N(N,q_new,dl,E,A,I,M,C,dti)
JE = zeros(2*N,2*N);
for i = 1:N-1
    if i == N-1
        JEs1 = hessEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),dl,E*A);
        JE(2*i-1:2*i+2,2*i-1:2*i+2) = JE(2*i-1:2*i+2,2*i-1:2*i+2)+JEs1;
        break
    end
    JEs1 = hessEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),dl,E*A);
    JEb = hessEb(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),q_new(2*i+3),q_new(2*i+4),0,dl,E*I);
    JE(2*i-1:2*i+2,2*i-1:2*i+2) = JE(2*i-1:2*i+2,2*i-1:2*i+2)+JEs1;
    JE(2*i-1:2*i+4,2*i-1:2*i+4) = JE(2*i-1:2*i+4,2*i-1:2*i+4)+JEb;
end
for k = 1:2*N
    for j = 1:2*N
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