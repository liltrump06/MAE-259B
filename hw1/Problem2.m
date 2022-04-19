clc;clear;close all;
N=21;
l = 0.1;
dl = l/(N-1);
R1 = dl/10;
Rmid = 0.025;

rou_m = 7000;
rou_f = 1000;
miu = 1000;
E = 1e9;
r0 = 0.001;
A = r0^2*pi;
I = pi*r0^4/4;
M = zeros(2*N,2*N);
C = zeros(2*N,2*N);
W = zeros(2*N,1);
g=9.81;
for i = 1:2*N
    if i == (N-1)+1 || i == (N-1)+2
        M(i,i) = 4/3*pi*Rmid^3*rou_m;
        C(i,i) = 6*pi*miu*Rmid;
    else
        M(i,i) = 4/3*pi*R1^3*rou_m;
        C(i,i) = 6*pi*miu*R1;
    end
    if i == (N-1)+2
        W(i) = -4/3*pi*Rmid^3*(rou_m-rou_f)*g;
        continue
    elseif rem(i,2) == 0
        W(i) = -4/3*pi*R1^3*(rou_m-rou_f)*g;
    end
end
%% initialize
dti = 0.01;
for i = 1:2*N
    if rem(i,2) == 0
        q_old(i) = 0;
    else
        q_old(i) = (i+1)/2*dl-dl;
    end
end
q_new = q_old +0.0005;
qd_old = zeros(2*N,1);
%% implicit
qlist = [];
qdlist = [];
tlist = [];
for i = 1:(50/dti-1)
    t= i*dti;
    err = 10;
    while err >1e-5
        f = ComputeF_N(N,q_new,q_old,qd_old,M,C,W,dl,E,A,I,dti);
        J = ComputeJ_N(N,q_new,dl,E,A,I,M,C,dti);
        deltaX = J \ f';
        q_new = q_new - deltaX';
        err = sum(abs(f));
    end
    qd_new = (q_new - q_old) / dti;
    qlist = [qlist;q_new];
    qdlist = [qdlist;qd_new];
    tlist = [tlist;t];
    q_old = q_new;
    qd_old = qd_new;
   
end
%%
FONT = 'Arial';
FONTSIZE = 10;
pWidth = 40; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 0 0 0]/255; % colors
figure
for k = 1:i
    plot(qlist(k,[1:2:2*N-1]),qlist(k,[2:2:2*N]),'ro-');
    xlim([0 0.1]);
    ylim([-0.3 0]);
    drawnow
end
%%
figure
plot(tlist,qlist(:,N+1))
legend('y position of middle node')
xlabel('time /s')
ylabel('y position /m')
title('vertical position of middle node vs time')
figure
plot(tlist,qdlist(:,N+1))
legend('y velocity of middle node')
xlabel('time /s')
ylabel('y velocity /m')
title('vertical velocity of middle node vs time')

