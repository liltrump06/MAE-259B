clc;clear;close all;
addpath("func\")
%% Phase1 
global N l dl dt time1 E A I free_index dampratio Rc m
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
        J = J + Jini+Jdamp;
        deltaX = J(free_index,free_index) \ F(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(F(free_index)));
    end
    qd_new = (q_new - q_old) / dt;
    qlist = [qlist,q_new];
    qrlist = zeros(size(qlist));
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    q_old = q_new;
    qd_old = qd_new;
end
%%
for o = 1: size(qlist,2)
    for k = 1:N
        qrlist(2*k-1:2*k,o) = Rc * qlist(2*k-1:2*k,o);
    end
end
save('phase1_rot.mat','q_new','q_old','qd_old','M','qrlist','dt','time1')
%%
qx = zeros(N,round(time1/dt));
qy = zeros(N,round(time1/dt));
for k =1:N
    qx(k,:) = qrlist(2*k-1,:);
    qy(k,:) = qrlist(2*k,:);
end
    
%v = VideoWriter('Output0.avi');
%open(v);
figure(1)
%set(gcf,'outerposition',get(0,'screensize'));
for k = 1:i
    if mod(k,5) == 0
        plot(qx(:,k),qy(:,k),'ro-');
        axis([0,0.2,-0.1,0.2])
        axis equal
    end
    
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    drawnow
end