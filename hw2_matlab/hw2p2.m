clc,clear,close all
N = 50;
q = zeros(4*N-1,1);
qd = zeros(4*N-1,1);
kappabarlist = zeros([N-2,2]);
dt = 0.01;
l = 0.2;
Rn = 0.02;
r0 = 0.001;
dtheta = l/Rn*1/(N-1);
rou = 1000;
E = 1e7;
G = E/3;
g = [0,0,-9.81];
J = pi*(r0^4)/2;
I = pi*(r0^4)/4;
A = pi*r0^2;
EA = E*A;
EI = E*I;
GJ = G*J;
M =eye(4*N-1);
M = M*(A*l*rou/(N-1));
for i = 1:N
    if i == N
        xk = [Rn* cos((i-1)*dtheta), Rn*sin((i-1)*dtheta),0];
        q(4*i-3:4*i-1) = xk';
        M(4*i-3:4*i-1,4*i-3:4*i-1) = diag(A*l*rou/(N)*ones(1,3));
    
    else
        xk = [Rn* cos((i-1)*dtheta), Rn*sin((i-1)*dtheta),0];
        q(4*i-3:4*i-1) = xk';
        q(4*i) = 0;
        M(4*i-3:4*i-1,4*i-3:4*i-1) = diag(A*l*rou/(N)*ones(1,3));
        M(4*i,4*i) = A*l*rou/(N-1)*r0^2/2;
    end
end
l_k = (norm(q(9:11)-q(5:7))+norm(q(5:7)-q(1:3)))/2;

%parameter prepare
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
m1 = zeros(N-1,3);
m2 = zeros(N-1,3);
mkref = zeros(N-2,1);
err = 100;
free_index = 8:4*N-1;
t1 = (q(5:7) - q(1:3))/norm((q(5:7) - q(1:3)));
u11 = [-t1(2),t1(1),0];
time = 5;
q_whole = zeros(4*N-1,time/dt);
qd_whole = zeros(4*N-1,time/dt);
kappanew = zeros(N-2,2);
for i = 1:N-2
    kappanew(i,:) = [0.2047929168960332,0];
end
for j = 1:time/dt
    err = 100;
    while err >1e-5
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
        [F,Jaco] = computeF_J(N,q_new,l_k,EA,m1,m2,GJ,EI,kappanew,mkref,g,M);
        F = F+ M*(((q_new-q_old)/dt-qd)/dt);
        Jini = M*(1/dt)*(1/dt);
        Jaco = Jaco + Jini;
        F_free = F(free_index);
        J_free = Jaco(free_index,free_index);
        Jinv=inv(J_free);
        deltaq = Jinv *F_free;
        q_new(free_index) = q_new(free_index) - deltaq ;
        err = sum(abs(F_free))
    end
    %figure(1)
    %for i = 1:N
    %    plot(i,q_new(4*i-1),'*')
    %    hold on 
    %end
    q_whole(:,j) = q_new;
    qd = (q_new-q_old)*(1/dt);
    q_old = q_new;
    a1_old = a1_new;
    a2_old = a2_new;
    tk_old = tk_new;
end
%%
plot(q_whole(4*N-1,:),'o')