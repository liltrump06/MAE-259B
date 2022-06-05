%% Phase2 prepare
clc;close all; clear
global N l dl dt E A I free_index stack stack_bound thumb dampratio Rc
Mat_prop;% load constants
dt = 0.0001; % set time step
free_index=1:2*N;
% read data from last phase
Phase1 = load('phase1_rot.mat');
q_new = Phase1.q_new;
q_old = Phase1.q_old;
qd_old = Phase1.qd_old; 
qr_new = zeros(size(q_new));
qrd_old = qr_new;
qr_old = qr_new;
M = Phase1.M;
time1 = Phase1.time1;
qd_new = qd_old;
e_dis = zeros(2*N,1);
qlist = [];
qdlist= [];
tlist = [];
stackx = zeros(50,1);
% perform the rotation to make coordinates in real configuration
for i = 1:N
    qr_new(2*i-1:2*i) = Rc * q_new(2*i-1:2*i);
    qr_old(2*i-1:2*i) = Rc * q_old(2*i-1:2*i);
    qrd_old(2*i-1:2*i) = Rc * qd_old(2*i-1:2*i);
end 
% compute the stack position using polyfit
[~,stack,thumb,~,stack_bound] = compute_stack(q_new);
% initialize constraint map
conmap = zeros(N,2); 
% initialize unit constraint direction vector
nlist = zeros(N,2);
% stack plot data
stackpl = zeros(5000,2);
stacky = linspace(q_new(2),q_new(2*N),5000);
for i = 1:5000
    stackx(i) = stacky(i)^2*stack(1)+stacky(i)*stack(2)+stack(3);
    stackpl(i,:) = Rc*[stackx(i);stacky(i)];
end

%% Phase2
time2 = 0.2;
fprintf('%s','enter phase2 \n')
nall = zeros(N,2,round(time2/dt+1));
for i = 1:(time2/dt+1)
    fprintf('time is %.4f s \n',time1 + i*dt)
    t= i*dt;
    err = 10;
    [z,M_modinv,nlist,interlist] = mass_mod_and_z_p2rot(qr_new,M,conmap,qrd_old);
%% predictor
    % First do a predict considering if the constraints are same as last
    % time step.
    while err >1e-5
        % only compute elastic force, gravity and its Jacobian
        [F,J] = ComputeF_J(N,qr_new,qr_old,qrd_old,M,e_dis);
        % use new Newton-Raphson equation to adapt to 
        % Mass modification method
        f = (qr_new-qr_old)/dt-qrd_old+dt*M_modinv*F-z;
        % compute reaction force
        Forceall = M*((qr_new-qr_old)/dt-qrd_old)/dt+F;
        % Jacobian in new method
        J = M_modinv*(dt)*J + eye(size(J))/dt + dampratio*M_modinv*eye(size(J));
        deltaX = J(free_index,free_index) \ f(free_index);
        qr_new(free_index) = qr_new(free_index) - deltaX;
        err = sum(abs(f(free_index)));
    end
    qrd_new =(qr_new-qr_old)/dt;
 %% corrector
 % use the predicted position and velocity to compute 
    % reset the calculation flag
    recal = 0;
    % change the constraint map according to the results of prediction
    [conmap,velop,recal] = judge_and_change_p2rot(qrd_new,qr_new,qr_old,conmap,Forceall);
    if recal == 1
        [z,M_modinv,nlist,interlist] = mass_mod_and_z_p2rot(qr_new,M,conmap,qrd_old);
        err=10;
        while err >1e-5
            [F,J] = ComputeF_J(N,qr_new,qr_old,qrd_old,M,e_dis);
            f = (qr_new-qr_old)/dt-qrd_old+dt*M_modinv*F-z;
            Forceall = M*((qr_new-qr_old)/dt-qrd_old)/dt+F;
            J = M_modinv*(dt)*J + eye(size(J))/dt+ dampratio*M_modinv*eye(size(J));
            deltaX = J(free_index,free_index) \ f(free_index);
            qr_new(free_index) = qr_new(free_index) - deltaX;
            err = sum(abs(f(free_index)));
        end
    
    end
    % update
    qlist = [qlist,qr_new];
    qrd_old = (qr_new-qr_old)/dt;
    qdlist = [qdlist,qrd_old];
    tlist = [tlist,t];
    nall(:,:,i) = nlist;
    qr_old = qr_new;
end
%% plot
qx = zeros(N,size(qlist,2));
qy = zeros(N,size(qlist,2));
for k = 1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end
figure(1)
set(gcf,'outerposition',get(0,'screensize'));
tic
for k = 1:i
    if mod(k,5) ==0
    plot(stackpl(:,1),stackpl(:,2),'-',qx(:,k),qy(:,k),'ro-');    
    drawnow
    end
end
toc
save('phase2_rot.mat','qr_new','qr_old','qrd_old','M','qlist','dt','time2')


