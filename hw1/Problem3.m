clc;clear;close all
N = 50;
l = 1;
dl = l/(N-1);
R = 0.013;
r = 0.011;
P = -2000;
E = 70e9;
I = pi/4*(R^4-r^4);
A = pi*(R^2-r^2);
rou_al = 2700;
m  = pi*(R^2-r^2)*l*rou_al/(N-1);
d = 0.75;
dti = 0.01;
for i = 1:2*N
    if rem(i,2) == 0
        q_old(i) = 0;
    else
        q_old(i) = (i+1)/2*dl-dl;
    end
end
%%
M = eye(N*2)*m;
Pv = zeros(2*N,1);
Pv(round((3/4)*2*N)+1) = P;
%%
free_index=3:(2*N-1);
fix_index = [1,2,2*N];
qlist = [];
qdlist= [];
tlist = [];
q_new = q_old;
qd_old = zeros(2*N,1);
for i = 1:(1/dti-1)
    t= i*dti;
    err = 10;
    while err >1e-5
        f = ComputeF_fix(N,q_new,q_old,qd_old,M,Pv,dl,E,A,I,dti);
        J = ComputeJ_fix(N,q_new,dl,E,A,I,M,dti);
        deltaX = J(free_index,free_index) \ f(free_index)';
        q_new(free_index) = q_new(free_index) - deltaX';
        err = sum(abs(f(free_index)));
    end
    q_new([1,2,2*N]) = 0;
    qd_new = (q_new - q_old) / dti;
    qlist = [qlist;q_new];
    qdlist = [qdlist;qd_new];
    tlist = [tlist;t];
    q_old = q_new;
    qd_old = qd_new;
end
%%
figure(2)
tic
for k = 1:i
    PL1 = plot(qlist(k,[1:2:(2*N-1)]),qlist(k,[2:2:2*N]),'ro-');
    axis([0 1 -0.4 0])
    drawnow
end
toc
%%
mlist = [];
for k = 1:size(qlist,1)
    m = min(qlist(k,:));
    mlist = [mlist;m];
end
figure
plot(tlist,mlist)
title('time vs ymin of all nodes')
xlabel('time /s')
ylabel('y position /m')
%%
c = min(l-d,d);
y_est = (P*c*(l^2-c^2)^1.5) / (9*sqrt(3)*E*I*l);

