function f = ComputeF_N(N,q_new,q_old,qd_old,M,C,W,dl,E,A,I,dti)
dE = zeros(2*N,1);
for i = 1:N-1
    if i == N-1
        dEs1 = gradEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),dl,E*A);
        dE(2*i-1:2*i+2) = dE(2*i-1:2*i+2)+dEs1;
        break
    end
    dEs1 = gradEs(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),dl,E*A);
    dEb = gradEb(q_new(2*i-1),q_new(2*i),q_new(2*i+1),q_new(2*i+2),q_new(2*i+3),q_new(2*i+4),0,dl,E*I);
    dE(2*i-1:2*i+2) = dE(2*i-1:2*i+2)+dEs1;
    dE(2*i-1:2*i+4) = dE(2*i-1:2*i+4)+dEb;
end
for m = 1:2*N 
        f(m) = M(m,m)/dti*((q_new(m)-q_old(m))/dti-qd_old(m))+dE(m)+C(m,m)*(q_new(m)-q_old(m))/dti-W(m);
end
end