%% Phase2 prepare
clc;close all; clear
global N l dl dt time1 E A I free_index stack polyn stack_bound thumb dampratio
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
q = zeros(2*N,1); 
time1 = 10;
free_index=1:2*N;
Phase1 = load('phase1.mat');
q_new = Phase1.q_new;
q_old = Phase1.q_old;
qd_old = Phase1.qd_old; 
M = Phase1.M;
qd_new = qd_old;
e_dis = zeros(2*N,1);
qlist = [];
qdlist= [];
tlist = [];
dampratio = 0.00001;
stackx = zeros(50,1);
%%
[~,stack,thumb,~,stack_bound] = compute_stack(q_new);% compute the stack position using polyfit
conmap = zeros(N,2);
nlist = zeros(N,2);
stacky = linspace(-0.04,0.2,5000);
for i = 1:5000
    stackx(i) = stacky(i)^2*stack(1)+stacky(i)*stack(2)+stack(3);
end
%% Phase2
time2 = 0.15;
fprintf('%s','enter phase2 \n')
nall = zeros(N,2,round(time2/dt+1));
%% predict
for i = 1:(time2/dt+1)
    fprintf('time is %.4f s \n',i*dt)

    t= i*dt;
    err = 10;
    [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z(q_new,M,conmap,qd_old);

    while err >1e-5
        %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
        f = (q_new-q_old)/dt-qd_old+dt*M_modinv*F-imposeacc;
        Forceall = M*((q_new-q_old)/dt-qd_old)/dt+F;
        J = M_modinv*(dt)*J + eye(size(J))/dt + dampratio*M_modinv*eye(size(J));
        deltaX = J(free_index,free_index) \ f(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(f(free_index)));
    end
 %% corrector
    recal = 0;
    [conmap,velop,recal] = judge_and_change(qd_new,q_new,q_old,conmap,Forceall);
    if recal == 1
        [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z(q_new,M,conmap,qd_old);

        err=10;
        %%
        for m = 1:N
            qx(m,:) = q_new(2*m-1);
            qy(m,:) = q_new(2*m);
        end

        %figure(3)
        
        %plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
        %axis([-0.0005,0.0005,0.0235,0.0245])

        %%
        while err >1e-5
            %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
            [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
            f = (q_new-q_old)/dt-qd_old+dt*M_modinv*F-imposeacc;
            Forceall = M*((q_new-q_old)/dt-qd_old)/dt+F;
            J = M_modinv*(dt)*J + eye(size(J))/dt+ dampratio*M_modinv*eye(size(J));
            deltaX = J(free_index,free_index) \ f(free_index);
            q_new(free_index) = q_new(free_index) - deltaX;
            err = sum(abs(f(free_index)));
        end
    end
    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    nall(:,:,i) = nlist;
    qd_old = (q_new-q_old)/dt;
    q_old = q_new;
    qx = zeros(N,1);
    qy = zeros(N,1);
    for m = 1:N
        qx(m,:) = q_new(2*m-1);
        qy(m,:) = q_new(2*m);
    end
    figure(1)
    plot(stackx,stacky,'-',qx,qy,'ro-')
    hold on
    quiver(qx,qy,nlist(:,1),nlist(:,2))
    hold off
    %figure(3)
    %plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
    %axis([-0.0005,0.0005,0.0235,0.0245])
end
%%
qx = zeros(N,size(qlist,2));
qy = zeros(N,size(qlist,2));
for k =1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end
    
%v = VideoWriter('Output0.avi');
%open(v);
figure(1)
%set(gcf,'outerposition',get(0,'screensize'));
tic
for k = 1:i
    if mod(k,10) ==0
    plot(stackx,stacky,'-',qx(:,k),qy(:,k),'ro-');
    hold on 
    quiver(qx(:,k),qy(:,k),nall(:,1,i),nall(:,2,i))
    hold off
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    drawnow
    end
end
toc


