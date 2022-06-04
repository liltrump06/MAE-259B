clc;clear;close all
%% Phase 3
addpath("func\")
global N l dl dt E A I free_index dampratio Rc ground M miu
Mat_prop;
dt = 0.001;
free_index=1:2*N;
%
s = load('Phase2_rot.mat');
q_new = s.qr_new;
q_old = s.qr_old;
qd_old = s.qrd_old;
M = s.M;
%
e_dis = zeros(2*N,1);
qlist = [];
qdlist= [];
tlist = [];
stackx = zeros(5000,1);
qd_new = qd_old;
Forceall = zeros(2*N,1);
%% constraint
ground = -0.3;
groundx = linspace(-0.3,1,1000);
groundy = ground.*ones(1,1000);
conmap = zeros(N,1);
%%
time2 = 0.6;
fprintf('%s','enter phase2 \n')
nall = zeros(N,2,round(time2/dt+1));
t = 0;
for i = 1:(time2/dt+1)
%     if sum(conmap) == 0
%         dt = 0.001;
%     else
%         dt = 0.0001;
%     end
    fprintf('time is %.4f s \n',t+dt)
    t=t+dt;
    err = 10;
    [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z_p3(q_new,M,conmap,qd_old);
%% predictor
    while err >1e-5
        %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
        [Ff,~] = compute_Ff_Jf(Forceall,conmap,J,M_modinv,qd_old);
        Forceall = (M*((q_new-q_old)/dt-qd_old)/dt+F);
        f = (q_new-q_old)/dt-qd_old+dt*M_modinv*(F-Ff)-imposeacc;        
        J = M_modinv*(dt)*J + eye(size(J))/dt + dampratio*M_modinv*eye(size(J));
        deltaX = J(free_index,free_index) \ f(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(f(free_index)));
    end
    qd_new =(q_new-q_old)/dt;
 %% corrector
    recal = 0;
    [conmap,recal] = judge_and_change_p3(qd_new,q_new,q_old,conmap,Forceall);
     if recal == 1
%         if sum(conmap) == 0
%             dt = 0.001;
%         else
%             dt = 0.0001;
%         end
        [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z_p3(q_new,M,conmap,qd_old);
        err=10;   
        while err >1e-5
            %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
            [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
            [Ff,~] = compute_Ff_Jf(Forceall,conmap,J,M_modinv,qd_old);
            Forceall = (M*((q_new-q_old)/dt-qd_old)/dt+F);
            f = (q_new-q_old)/dt-qd_old+dt*M_modinv*(F-Ff)-imposeacc;
            J = M_modinv*(dt)*J + eye(size(J))/dt + dampratio*M_modinv*eye(size(J));
            deltaX = J(free_index,free_index) \ f(free_index);
            err = sum(abs(f(free_index)));
            q_new(free_index) = q_new(free_index) - deltaX;
        end
    
    end
    qlist = [qlist,q_new];
    qd_old = (q_new-q_old)/dt;
    qdlist = [qdlist,qd_old];
    nall(:,:,i) = nlist;
    q_old = q_new;
    qx = zeros(N,1);
    qy = zeros(N,1);
    for m = 1:N
        qx(m,:) = q_new(2*m-1);
        qy(m,:) = q_new(2*m);
    end
    figure(1)
    set(gcf,'outerposition',get(0,'screensize'));
    plot(groundx,groundy,'-',qx,qy,'ro-')
    hold on
    quiver(qx,qy,nlist(:,1),nlist(:,2))
    hold off
    axis equal
%     figure(3)
%     plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
%     axis([-0.0005,0.0005,0.0235,0.0245])
end
%%
qx = zeros(N,size(qlist,2));
qy = zeros(N,size(qlist,2));
for k = 1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end
    
%v = VideoWriter('Output0.avi');
%open(v);
figure(1)
%set(gcf,'outerposition',get(0,'screensize'));
tic
for k = 1:i
    if mod(k,1) ==0
        plot(groundx,groundy,'-',qx(:,k),qy(:,k),'ro-');
        axis([0.2,0.8,-0.6,0.2])
        drawnow
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    end
end
toc
save('phase3.mat','qlist','dt','time2')


