function f = ComputeF(q_new,q_old,qd_old,M,C,W,l,E,A,I,dti)

dEs1 = gradEs(q_new(1),q_new(2),q_new(3),q_new(4),l/2,E*A);
dEs2 = gradEs(q_new(3),q_new(4),q_new(5),q_new(6),l/2,E*A);
dEb = gradEb(q_new(1),q_new(2),q_new(3),q_new(4),q_new(5),q_new(6),0,l/2,E*I);
for m = 1:6 
    if m == 1 ||m == 2
        f(m) = M(m,m)/dti*((q_new(m)-q_old(m))/dti-qd_old(m))+dEs1(m)+dEb(m)+C(m,m)*(q_new(m)-q_old(m))/dti-W(m,m);
    elseif m == 3 || m == 4
        f(m) = M(m,m)/dti*((q_new(m)-q_old(m))/dti-qd_old(m))+dEs1(m)+dEb(m)+dEs2(m-2)+C(m,m)*(q_new(m)-q_old(m))/dti-W(m,m);
    else
        f(m) = M(m,m)/dti*((q_new(m)-q_old(m))/dti-qd_old(m))+dEb(m)+dEs2(m-2)+C(m,m)*(q_new(m)-q_old(m))/dti-W(m,m);
    end
end
end