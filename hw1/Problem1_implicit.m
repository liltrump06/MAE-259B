clc;clear;
R1 = 0.005;
R2 = 0.025;
R3 = 0.005;
rou_m = 7000;
rou_f = 1000;
miu = 1000;
l = 0.1;
E = 1e9;
r0 = 0.001;
A = r0^2*pi;
I = pi*r0^4/4;
M = zeros(6,6);
M(1,1) = 4/3*pi*R1^3*rou_m;
M(2,2) = 4/3*pi*R1^3*rou_m;
M(3,3) = 4/3*pi*R2^3*rou_m;
M(4,4) = 4/3*pi*R2^3*rou_m;
M(5,5) = 4/3*pi*R3^3*rou_m;
M(6,6) = 4/3*pi*R3^3*rou_m;
C = zeros(6,6);
C(1,1) = 6*pi*miu*R1;
C(2,2) = 6*pi*miu*R1;
C(3,3) = 6*pi*miu*R2;
C(4,4) = 6*pi*miu*R2;
C(5,5) = 6*pi*miu*R3;
C(6,6) = 6*pi*miu*R3;
W = zeros(6,1);
g=9.81;
W(2,2) = -4/3*pi*R1^3*(rou_m-rou_f)*g;
W(4,4) = -4/3*pi*R2^3*(rou_m-rou_f)*g;
W(6,6) = -4/3*pi*R3^3*(rou_m-rou_f)*g;

%% implicit
dti = 0.01;
q_old = [0,0,0.05,0,0.1,0]';
qd_old = zeros(6,1);
q_new = [0.0005,-0.0005,0.05,-0.003,0.0995,-0.0005]';
err = 1;
qlist = [];
qdlist= [];
for i = 1:(10/dti-1)
    t= i*dti;
    err = 10;
    while err >1e-5
        f = ComputeF(q_new,q_old,qd_old,M,C,W,l,E,A,I,dti);
        J = ComputeJ(q_new,l,E,A,I,M,C,dti);
        deltaX = J \ f';
        q_new = q_new - deltaX;
        err = sum(abs(f));
    end
    qd_new = (q_new - q_old) / dti;
    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    q_old = q_new;
    qd_old = qd_new;
end
%%
FONT = 'Arial';
FONTSIZE = 10;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 0 0 0]/255; % colors
figure
for k = 1:i
    plot(qlist([1,3,5],k),qlist([2,4,6],k),'ro-');
    xlim([0 0.1]);
    ylim([-0.1 0]);
    drawnow
end



