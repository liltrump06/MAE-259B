%% Phase2 prepare
addpath("func\")
clc;close all; clear
global N l dl dt time1 E A I free_index stack polyn stack_bound thumb dampratio 
dt = 0.0001;
Mat_prop;
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
stackx = zeros(50,1);
%% 
[~,stack,thumb,~,stack_bound] = compute_stack(q_new);% compute the stack position using polyfit
conmap = zeros(N,2);
nlist = zeros(N,2);
stacky = linspace(q_new(2),q_new(2*N),5000);
for i = 1:5000
    stackx(i) = stacky(i)^2*stack(1)+stacky(i)*stack(2)+stack(3);
end
thumbx = linspace(-0.003,0.02,500);
thumby = thumb*ones(1,500);
%% Phase2
time2 = 0.08;
fprintf('%s','enter phase2 \n')
nall = zeros(N,2,round(time2/dt+1));
v = VideoWriter('Phase2.avi');
open(v);
figure(1)
for i = 1:(time2/dt+1)
    fprintf('time is %.4f s \n',i*dt)
    t= i*dt;
    err = 10;
    [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z(q_new,M,conmap,qd_old);
%% predict
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
    qd_new = (q_new-q_old)/dt;
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
%         figure(1)
%         plot(stackx,stacky,'-',qx,qy,'ro-')
%         hold on
%         quiver(qx,qy,nlist(:,1),nlist(:,2))
%         hold off
%         axis([-0.1,0.3,-0.05,0.3])
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
    qd_new = qd_old;
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
    set(gcf,'outerposition',get(0,'screensize'));
    plot(stackx,stacky,'b-',thumbx,thumby,'b-','LineWidth',2)
    hold on
    plot(qx,qy,'ro-')
    quiver(qx,qy,nlist(:,1),nlist(:,2))
    hold off
    axis([-0.1,0.3,-0.05,0.3])
    frame = getframe(gcf);
    writeVideo(v,frame);
    %figure(3)
    %plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
    %axis([-0.0005,0.0005,0.0235,0.0245])
end
save('phase2.mat','q_new','q_old','qd_old','M','qlist','dt','time2')
% %%
% qx = zeros(N,size(qlist,2));
% qy = zeros(N,size(qlist,2));
% for k =1:N
%     qx(k,:) = qlist(2*k-1,:);
%     qy(k,:) = qlist(2*k,:);
% end
% 
% %v = VideoWriter('Output0.avi');
% %open(v);
% figure(1)
% %set(gcf,'outerposition',get(0,'screensize'));
% tic
% for k = 1:i
%     if mod(k,5) ==0
%         plot(stackx,stacky,'b-',thumbx,thumby,'b-',qx(:,k),qy(:,k),'ro-');
%         drawnow
%         %frame = getframe(gcf);
%         %writeVideo(v,frame);
%     end
% end
% toc


