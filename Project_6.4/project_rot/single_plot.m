clc;clear;close all
s1 = load('phase1_rot.mat');
s2 = load('phase2_rot.mat');
s3 = load('phase3.mat');
%%
N = 50;
q1 = s1.qrlist;
q2 = s2.qlist;
q3 = s3.qlist;
%% resample
dt2 = s2.dt;
dt1 = s1.dt;
q2n = [];
for i = 1:size(q2,2)
    if mod(i,dt1/dt2) == 0
        q2n = [q2n,q2(:,i)];
    end
end
qall = [q1,q2n,q3];
time = size(qall,2);
qx = zeros(N,time);
qy = zeros(N,time);
%%
for i = 1:N
    qx(i,:) = qall(2*i-1,:);
    qy(i,:) = qall(2*i,:);
end
qz = zeros(size(qy));
qx50_1 = qx(50,:);
qy50_1 = qy(50,:);
qz50_1 = 0.05*ones(1,time);
qx50_2 = qx(50,:);
qy50_2 = qy(50,:);
qz50_2 = 0.1*ones(1,time);
qx50_3 = qx(50,:);
qy50_3 = qy(50,:);
qz50_3 = 0.15*ones(1,time);
qx1_1 = qx(1,:);
qy1_1 = qy(1,:);
qz1_1 = 0.15*ones(1,time);
qx1_2 = qx(1,:);
qy1_2 = qy(1,:);
qz1_2 = 0.1*ones(1,time);
qx1_3 = qx(1,:);
qy1_3 = qy(1,:);
qz1_3 = 0.05*ones(1,time);
qxall = [qx;qx50_1;qx50_2;qx50_3;flip(qx,1);qx1_1;qx1_2;qx1_3;qx(1,:)];
qyall = [qy;qy50_1;qy50_2;qy50_3;flip(qy,1);qy1_1;qy1_2;qy1_3;qy(1,:)];
qzall = [qz;qz50_1;qz50_2;qz50_3;flip(qz,1)+0.15;qz1_1;qz1_2;qz1_3;qz(1,:)];
insertx1 = qx(:,size(q1,2));
inserty1 = qy(:,size(q1,2));
insertz1 = qz(:,size(q1,2));
qxall()
%%
figure
for i = 1:time
    if mod(i,3) == 0
        plot3(qx(:,i),qz(:,i),qy(:,i),'ro-','LineWidth',0.5)
        hold on
        plot3(qx(:,i),qz(:,i)+0.05,qy(:,i),'ro-','LineWidth',0.5)
        plot3(qx(:,i),qz(:,i)+0.1,qy(:,i),'ro-','LineWidth',0.5)
        plot3(qx(:,i),qz(:,i)+0.15,qy(:,i),'ro-','LineWidth',0.5)
        hold off
        axis equal
        view([30,10])
        grid on
        drawnow
    end
end
%%
figure
for i = 1:time
    if mod(i,3) == 0
        plot3(qxall(:,i),qzall(:,i),qyall(:,i),'r-','LineWidth',2)
        hold on
        plot3(qxall(:,i)-0.02,qzall(:,i),qyall(:,i),'r-','LineWidth',2)
        plot3(qxall(:,i)-0.04,qzall(:,i),qyall(:,i),'r-','LineWidth',2)
        plot3(qxall(:,i)-0.06,qzall(:,i),qyall(:,i),'r-','LineWidth',2)
        hold off
        axis([-0.1,0.6,-0.3,0.4,-0.35,0.35])
        view([30,10])
        grid on
        drawnow
    end
end
%%
v = VideoWriter('single_plot.avi');
open(v);
figure(1)
set(gcf,'outerposition',get(0,'screensize'));
for i = 1:time
    if mod(i,2) == 0
        plot3(qxall(:,i),qzall(:,i),qyall(:,i),'r-','LineWidth',2)
        axis([qx(1,i)-0.15,qx(1,i)+0.15,qz(1,i)-0.05,qz(1,i)+0.25,qy(1,i)-0.05,qy(1,i)+0.25])
        view([30,10])
        grid on
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
