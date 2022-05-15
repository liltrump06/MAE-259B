tic
clc,clear,close all
global N l

N = 50; % define the number of nodes
q = zeros(4*N-1,1);% initialize the DOF vector
qd = zeros(4*N-1,1);
dt = 0.001;% set the time interval
%% basic property value of the rod
l = 0.14605; % length
r0 =  2.9074e-04/2; % radius of the rod
rou = 436; % density kg/m3
E = 1e9; % Young's Modulus
G = E/3; % shear modulus ,usually equal to E/3 when the material is elastic.
g = [0,0,-9.81]; % gravity acceleration
J = pi*(r0^4)/2; 
I = pi*(r0^4)/4; 
A = pi*r0^2; % area of the cross section
EA = E*A;
EI = E*I;
GJ = G*J;
M = eye(4*N-1);% initialize the mass matrix
M = M*(A*l*rou/(N-1));% distribute the mass to each node
%% initialize q at t=0
cood_ini = linspace(0,l,N);
for i = 1:N
    if i == N % there are only N-1 theta(s)
        xk = [0,0,cood_ini(i)];
        q(4*i-3:4*i-1) = xk';
        M(4*i-3:4*i-1,4*i-3:4*i-1) = diag(A*l*rou/(N)*ones(1,3)); % distribute the mass to each node
    
    else
        xk = [0,0,cood_ini(i)];
        q(4*i-3:4*i-1) = xk';
        q(4*i) = 0;
        M(4*i-3:4*i-1,4*i-3:4*i-1) = diag(A*l*rou/(N)*ones(1,3));
        M(4*i,4*i) = A*l*rou/(N-1)*r0^2/2; %calculate the inertia
    end
end
l_k = (norm(q(9:11)-q(5:7))+norm(q(5:7)-q(1:3)))/2; % vonoroi length
M_mod = M;
%% parameter initialize
u1 = zeros(N-1,3);
u2 = zeros(N-1,3);
tk_old = zeros(N-1,3);
tk_new = zeros(N-1,3);
a1_old = zeros(N-1,3);
a1_new = zeros(N-1,3);
a2_old = zeros(N-1,3);
a2_new = zeros(N-1,3);
q_new = q;
q_old = q;
q_realold = zeros(size(q));
m1 = zeros(N-1,3);
m2 = zeros(N-1,3);
mkref = zeros(N-2,1);
err = 100;
free_index = 3:4*(N-1)-1;
t1 = (q(5:7) - q(1:3))/norm((q(5:7) - q(1:3)));
u11 = [1,0,0];
time = 5;
q_whole = zeros(4*N-1,time/dt);
qd_whole = zeros(4*N-1,time/dt);
kappanew = zeros(N-2,2); % initialize kappabar vector in undeformed status
e_dis = zeros(4*N-1,1);
%% get reference frame and material frame at t=0 to calculate undeformed kappa
for m = 1:N-1
    tk1 =(q_new(4*m+1:4*m+3) - q_new(4*m-3:4*m-1))/norm((q_new(4*m+1:4*m+3) - q_new(4*m-3:4*m-1)));
    if m == 1
       u1(m,:) = u11;
    else
       u1(m,:) = parallel_transport(u1(m-1,:)',tk_new(m-1,:)',tk1);
    end
    a1_new(m,:) = u1(m,:);      
    a2_new(m,:)= cross(tk1,a1_new(m,:));
    u2(m,:)= cross(tk1,u1(m,:));
    tk_new(m,:)= tk1;
    m1(m,:)= cos(q_new(4*m))*a1_new(m,:)+sin(q_new(4*m))*a2_new(m,:);
    m2(m,:)= cos(q_new(4*m))*a2_new(m,:)-sin(q_new(4*m))*a1_new(m,:);
end
for i = 1:N-2
    ka = computekappa(q_new(4*i-3:4*i-1),q_new(4*i+1:4*i+3),q_new(4*i+5:4*i+7), ...
        m1(i,:),m2(i,:),m1(i+1,:),m2(i+1,:));
    kappanew(i,:) = ka;
end
% create the external_force distribution
%e_dis = compute_external_force_distribution(jlist,N);
%% LOOP
for j = 1:(time/dt)% time loop
    err = 100;
    %free_index change 
    if j == time/dt/5+1 % 1001 time point
        free_index = 1:4*N-1;% release the rod when the external force = 0
    end
    while err >1e-6 % error loop
        for m = 1:N-1
            tk1 =(q_new(4*m+1:4*m+3) - q_new(4*m-3:4*m-1))/norm((q_new(4*m+1:4*m+3) - q_new(4*m-3:4*m-1)));
            if j == 1
                if m == 1
                    u1(m,:) = u11;
                else
                    u1(m,:) = parallel_transport(u1(m-1,:)',tk_new(m-1,:)',tk1);
                end
                a1_new(m,:) = u1(m,:);                
            else
                if m == 1
                    u1(m,:)= u11;
                else
                    u1(m,:)= parallel_transport(u1(m-1,:)',tk_new(m-1,:)',tk1);
                end
                a1_new(m,:)= parallel_transport(a1_old(m,:)',tk_old(m,:)',tk1);
            end
            a2_new(m,:)= cross(tk1,a1_new(m,:));
            u2(m,:)= cross(tk1,u1(m,:));
            tk_new(m,:)= tk1;
            m1(m,:)= cos(q_new(4*m))*a1_new(m,:)+sin(q_new(4*m))*a2_new(m,:);
            m2(m,:)= cos(q_new(4*m))*a2_new(m,:)-sin(q_new(4*m))*a1_new(m,:);
        end
        
        for k = 1:N-2
            mkref(k) = signedAngle(parallel_transport(a1_new(k,:)',tk_new(k,:)', ...
                tk_new(k+1,:)'),a1_new(k+1,:)',tk_new(k,:)');
        end
        %e_dis = compute_external_force_distribution(j,N,q_new,q_old,dt,time);
        [F,Jaco] = computeF_J(N,q_new,l_k,EA,m1,m2,GJ,EI,kappanew,mkref,g,M,e_dis);% compute F and J
        F = F+ M_mod*(((q_new-q_old)/dt-qd)/dt);
        Jini = M_mod*(1/dt)*(1/dt);
        Jaco = Jaco + Jini;
        F_free = F(free_index);
        J_free = Jaco(free_index,free_index);
        Jinv=inv(J_free);
        deltaq = Jinv *F_free;
        q_new(free_index) = q_new(free_index) - deltaq ;
        err = sum(abs(F_free));
    end
    %figure(1)
    %for i = 1:N
    %    plot(i,q_new(4*i-1),'*')
    %    hold on 
    %end
    % free_index change due to contact force
    % compute the constraint
    if j == time/dt/5
        [~,stack_con,thumb_con,polynum] = compute_stack(q_new);
    end
    if j > time/dt/5
        %contact_standard()
    end
    q_whole(:,j) = q_new;% store the q data at every time point.
    qd = (q_new-q_old)*(1/dt);
    e_dis = compute_external_force_distribution(j+1,q_new,q_old,qd,dt,time);
    q_old = q_new;
    a1_old = a1_new;
    a2_old = a2_new;
    tk_old = tk_new;
    j
end
q_whole = [q,q_whole];% cat qwhole with the initial state q.
toc
%%
figure(1)
plot(1:(time/dt+1),q_whole(3,:),'o')
%%
save('project.mat','q_whole','-mat')
timep = size(q_whole,2);
qplot = zeros(N,3,timep);
for i = 1:N
    for t = 1:timep
        qplot(i,:,t) = q_whole(4*i-3:4*i-1,t);% reshape the matrix
    end
end
qcard = qplot;
q_para = qplot;
for i = 1:5 % this loop create another five rod to simulate the card
    q_para(:,2,:) = qplot(:,2,:)+0.01*i;
    qcard = cat(1,qcard,q_para);
end
%%

%v = VideoWriter('Output2.avi');
%open(v);
figure(2)
%set(gcf,'outerposition',get(0,'screensize'));

for t = 1:size(q_whole,2)
    if mod(t,10) ==0
        scatter3(qcard(:,1,t),qcard(:,2,t),qcard(:,3,t),'o')% simulate the process
        axis([-0.1,0.3,-0.1,0.1,-0.2,0.2])
        view([10,10])
        %frame = getframe(gcf);
        %writeVideo(v,frame);
        drawnow
    end
end
