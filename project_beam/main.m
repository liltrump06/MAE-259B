clc;clear;close all;
global N l dl dt time1 E A I free_index stack polyn
N = 50;
l = 0.14605; % length
dl = l/(N-1);
r0 =  2.9074e-04/2; % radius of the rod
rou = 436; % density kg/m3
E = 1e9; % Young's Modulus
I = pi/4*(r0^4);
A = pi*r0^2;
rou_al = 436;
m  = pi*(r0^2)*l*rou_al/(N-1);
dt = 0.001;
q=zeros(2*N,1);
cood_ini = linspace(0,l,N);
for i = 1:N
    xk = [0,cood_ini(i)];
    q(2*i-1:2*i) = xk';
end
M = eye(N*2)*m;
time1 = 1.5;
free_index=2:2*(N-1);

qlist = [];
qdlist= [];
tlist = [];
q_new = q;
q_old = q;
qd_old = zeros(2*N,1);
%% Phase1
for i = 1:(time1/dt)
    fprintf('time is %.4f s \n',i*dt)
    if i == round(time1/dt/5)+1
        free_index = 1:2*N;
    end
    t= i*dt;
    err = 10;
    while err >1e-5
        M_mod = mass_mod(M,conmap,)
        e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,M,e_dis,dl,E,A,I);
        F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
        Jini = M*(1/dt)*(1/dt);
        J = J + Jini;
        deltaX = J(free_index,free_index) \ F(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(F(free_index)));
    end
    qd_new = (q_new - q_old) / dt;
    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    q_old = q_new;
    qd_old = qd_new;
end
%% Phase2 prepare
[~,stack,thumb,polyn] = compute_stack(q_new);% compute the stack position using polyfit
conmap = zeros(N,1);
nlist = zeros(N,2);
%% Phase2
time2 = 5;
fprintf('%s','enter phase2')
for i = 1:(time2/dt)
    fprintf('time is %.4f s \n',i*dt+time1)
    if i == round(time2/dt/5)+1
        free_index = 1:2*N;
    end
    t= i*dt+time1;
    err = 10;
    while err >1e-5


        e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,M,e_dis,dl,E,A,I);
        F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
        Jini = M*(1/dt)*(1/dt);
        J = J + Jini;
        deltaX = J(free_index,free_index) \ F(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(F(free_index)));
    end
    qd_new = (q_new - q_old) / dt;
    recal = 0;
    [conmap,velop,recal] = judge_and_change(qd_new,q_new,q_old,conmap,F);
    if recal == 1

        while err >1e-5


            e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
            [F,J] = ComputeF_J(N,q_new,M,e_dis,dl,E,A,I);
            F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
            Jini = M*(1/dt)*(1/dt);
            J = J + Jini;
            deltaX = J(free_index,free_index) \ F(free_index);
            q_new(free_index) = q_new(free_index) - deltaX;
            err = sum(abs(F(free_index)));
        end





    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    q_old = q_new;
    qd_old = qd_new;
end
%%
qx = zeros(N,time/dt);
qy = zeros(N,time/dt);
for k =1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end
    
figure
for k = 1:i
    plot(qx(:,k),qy(:,k),'ro-');
    axis([-0.1,0.1,-0.1,0.2])
    drawnow
end


