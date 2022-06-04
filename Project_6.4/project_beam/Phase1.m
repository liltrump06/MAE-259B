clc;clear;close all;
%%
addpath("func\")
global N l dl dt time1 E A I free_index dampratio 
dt = 0.001;
Mat_prop;
cood_ini = linspace(0,l,N);
for i = 1:N
    xk = [0,cood_ini(i)];
    q(2*i-1:2*i) = xk';
end
M = eye(N*2)*m;
time1 = 0.7;
free_index=2:2*(N-1);
qlist = [];
qdlist= [];
tlist = [];
elist = [];
q_new = q;
q_old = q;
qd_old = zeros(2*N,1);
%% Phase1
for i = 1:(time1/dt+1)
    fprintf('time is %.4f s \n',i*dt)
    t= i*dt;
    err = 10;
    while err >1e-5
        e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
        F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
        Jini = M*(1/dt)*(1/dt);
        Jdamp = dampratio/dt*eye(size(J));
        J = J + Jini + Jdamp;
        deltaX = J(free_index,free_index) \ F(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(F(free_index)));
    end
    elist = [elist,e_dis];
    qd_new = (q_new - q_old) / dt;
    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    q_old = q_new;
    qd_old = qd_new;
end
%
save('phase1.mat','q_new','q_old','qd_old','M','dt','time1','qlist')
%%
qx = zeros(N,round(time1/dt));
qy = zeros(N,round(time1/dt));
vx = zeros(N,round(time1/dt));
vy = zeros(N,round(time1/dt));
for k =1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
    vx(k,:) = elist(2*k-1,:);
    vy(k,:) = elist(2*k,:);
end
    
v = VideoWriter('phase1_with_PD.avi');
open(v);
figure(1)
set(gcf,'outerposition',get(0,'screensize'));
for k = 1:i
    if mod(k,3) == 1
        plot(qx(:,k),qy(:,k),'ro-');
        hold on
        quiver(qx(:,k),qy(:,k),vx(:,k),vy(:,k))
        axis([-0.1,0.1,-0.05,0.15])

        frame = getframe(gcf);
        writeVideo(v,frame);
        drawnow
        hold off
    end
end