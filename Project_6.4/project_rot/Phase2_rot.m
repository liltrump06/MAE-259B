%% Phase2 prepare
clc;close all; clear
global N l dl dt E A I free_index stack stack_bound thumb dampratio Rc
Mat_prop;
dt = 0.0001; 
free_index=1:2*N;
Phase1 = load('phase1_rot.mat');
q_new = Phase1.q_new;
q_old = Phase1.q_old;
qd_old = Phase1.qd_old; 
qr_new = zeros(size(q_new));
qrd_old = qr_new;
qr_old = qr_new;
M = Phase1.M;
qd_new = qd_old;
e_dis = zeros(2*N,1);
qlist = [];
qdlist= [];
tlist = [];
stackx = zeros(50,1);
for i = 1:N
    qr_new(2*i-1:2*i) = Rc * q_new(2*i-1:2*i);
    qr_old(2*i-1:2*i) = Rc * q_old(2*i-1:2*i);
    qrd_old(2*i-1:2*i) = Rc * qd_old(2*i-1:2*i);
end
%% 
[~,stack,thumb,~,stack_bound] = compute_stack(q_new);% compute the stack position using polyfit
conmap = zeros(N,2);
nlist = zeros(N,2);
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
    fprintf('time is %.4f s \n',i*dt)
    t= i*dt;
    err = 10;
    [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z_p2rot(qr_new,M,conmap,qrd_old);
%% predictor
    while err >1e-5
        %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,qr_new,qr_old,M,e_dis,dl,E,A,I);
        f = (qr_new-qr_old)/dt-qrd_old+dt*M_modinv*F-imposeacc;
        Forceall = M*((qr_new-qr_old)/dt-qrd_old)/dt+F;
        J = M_modinv*(dt)*J + eye(size(J))/dt + dampratio*M_modinv*eye(size(J));
        deltaX = J(free_index,free_index) \ f(free_index);
        qr_new(free_index) = qr_new(free_index) - deltaX;
        err = sum(abs(f(free_index)));
    end
    qrd_new =(qr_new-qr_old)/dt;
 %% corrector
    recal = 0;
    [conmap,velop,recal] = judge_and_change_p2rot(qrd_new,qr_new,qr_old,conmap,Forceall);
    if recal == 1
        [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z_p2rot(qr_new,M,conmap,qrd_old);

        err=10;
        %%
        for m = 1:N
            qx(m,:) = qr_new(2*m-1);
            qy(m,:) = qr_new(2*m);
        end
        figure(1)
        set(gcf,'outerposition',get(0,'screensize'));
        plot(stackpl(:,1),stackpl(:,2),'-',qx,qy,'ro-')
        hold on
        quiver(qx,qy,nlist(:,1),nlist(:,2))
        hold off
        axis equal
        %figure(3)
        
        %plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
        %axis([-0.0005,0.0005,0.0235,0.0245])

        %%
        while err >1e-5
            %e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
            [F,J] = ComputeF_J(N,qr_new,qr_old,M,e_dis,dl,E,A,I);
            f = (qr_new-qr_old)/dt-qrd_old+dt*M_modinv*F-imposeacc;
            Forceall = M*((qr_new-qr_old)/dt-qrd_old)/dt+F;
            J = M_modinv*(dt)*J + eye(size(J))/dt+ dampratio*M_modinv*eye(size(J));
            deltaX = J(free_index,free_index) \ f(free_index);
            qr_new(free_index) = qr_new(free_index) - deltaX;
            err = sum(abs(f(free_index)));
        end
    
    end
    qlist = [qlist,qr_new];
    qrd_old = (qr_new-qr_old)/dt;
    qdlist = [qdlist,qrd_old];
    tlist = [tlist,t];
    nall(:,:,i) = nlist;
    qr_old = qr_new;
    qx = zeros(N,1);
    qy = zeros(N,1);
    for m = 1:N
        qx(m,:) = qr_new(2*m-1);
        qy(m,:) = qr_new(2*m);
    end
    figure(1)
    set(gcf,'outerposition',get(0,'screensize'));
    plot(stackpl(:,1),stackpl(:,2),'-',qx,qy,'ro-')
    hold on
    quiver(qx,qy,nlist(:,1),nlist(:,2))
    hold off
    axis equal
    %figure(3)
    %plot(stackx,stacky,'-',interlist(:,1),interlist(:,2),'ro',qx,qy,'bo')
    %axis([-0.0005,0.0005,0.0235,0.0245])
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
    if mod(k,5) ==0
    plot(stackpl(:,1),stackpl(:,2),'-',qx(:,k),qy(:,k),'ro-');
    drawnow
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    end
end
toc
save('phase2_rot.mat','qr_new','qr_old','qrd_old','M','qlist','dt','time2')


