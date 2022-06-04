clc;clear;close all
addpath('func')
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
qinsert = q1(:,size(q1,2));
qin_list = [];
for m = 1:25
    qin_list = [qin_list,qinsert];
end

%%
 qall1 = [q1,q2n,q3];
 qall2 = [q1,qin_list,q2n,q3];
 qall3 = [q1,qin_list,qin_list,q2n,q3];
 qall4 = [q1,qin_list,qin_list,qin_list,q2n,q3];
% qall5 = [q1,qin_list,q2n,q3];
% qall6 = [q1,qin_list,qin_list,q2n,q3];
% qall7 = [q1,qin_list,qin_list,qin_list,q2n,q3];
% [qxall1,qyall1,qzall1] = make_card(qall1,N);
% [qxall2,qyall2,qzall2] = make_card(qall2,N);
% [qxall3,qyall3,qzall3] = make_card(qall3,N);
% [qxall4,qyall4,qzall4] = make_card(qall4,N);
% [qxall5,qyall5,qzall5] = make_card(qall5,N);
% [qxall6,qyall6,qzall6] = make_card(qall6,N);
% [qxall7,qyall7,qzall7] = make_card(qall7,N);
% figure
time = 1501;
% for i = 1:time
%     if mod(i,3) == 0
%         plot3(qxall1(:,i),qzall1(:,i),qyall1(:,i),'b-','LineWidth',2)
%         hold on
%         plot3(qxall2(:,i)-0.001,qzall2(:,i),qyall2(:,i),'b-','LineWidth',2)
%         plot3(qxall3(:,i)-0.002,qzall3(:,i),qyall3(:,i),'b-','LineWidth',2)
%         plot3(qxall4(:,i)-0.003,qzall4(:,i),qyall4(:,i),'b-','LineWidth',2)
%         plot3(qxall2(:,i)-0.001,qzall2(:,i),qyall2(:,i),'b-','LineWidth',2)
%         plot3(qxall3(:,i)-0.002,qzall3(:,i),qyall3(:,i),'b-','LineWidth',2)
%         plot3(qxall4(:,i)-0.003,qzall4(:,i),qyall4(:,i),'b-','LineWidth',2)
%         hold off
%         axis([-0.1,0.6,-0.3,0.4,-0.35,0.35])
%         view([30,10])
%         grid on
%         drawnow
%     end
% end
% %%
% figure
% for i = 1:time
%     if mod(i,3) == 0
%         plot3(qxall1(:,i),qzall1(:,i),qyall1(:,i),'r-','LineWidth',2)
%         axis([qxall1(1,i)-0.15,qxall1(1,i)+0.15,qzall1(1,i)-0.05,qzall1(1,i)+0.25,qyall1(1,i)-0.05,qyall1(1,i)+0.25])
%         view([30,10])
%         grid on
%         drawnow
%     end
% end
%%
original1 = imread('Ace_of_spades.png');
original2 = imread('2_of_diamonds.png');
original3 = imread('3_of_clubs.png');
original4 = imread('4_of_hearts.png');
original5 = imread('wood.jpg');
[qx1,qy1,qz1] = make_card_surf(qall1,N);
[qx2,qy2,qz2] = make_card_surf(qall2,N);
[qx3,qy3,qz3] = make_card_surf(qall3,N);
[qx4,qy4,qz4] = make_card_surf(qall4,N);
figure(1)
%%
for i = 1:time
    if mod(i,3) == 0
        qxsurf = [];
        qzsurf = [];
        for l = 1:N
            qxsurf = [qxsurf,qx1(:,i)];            
            qzsurf = [qzsurf;linspace(0,0,50)];
        end
        [mg1,mg2] = meshgrid(linspace(0,0.2,50),qy1(:,i));
        m = surf(mg1,mg2,qxsurf);
        axis([qz1(1,i)-0.05,qz1(1,i)+0.25,qy1(1,i)-0.05,qy1(1,i)+0.25,qx1(1,i)-0.15,qx1(1,i)+0.15,])
        view([80,170])
        grid on
        camroll(270)
        h = findobj('Type', 'surface');
        set(h, 'CData', original1, 'FaceColor', 'texturemap')
        set(m, 'edgecolor','none')
        set(gcf,'outerposition',get(0,'screensize'));
        drawnow
    end  
end

%%
v = VideoWriter('multi_plot.avi');
open(v);
figure(2)
set(gcf,'outerposition',get(0,'screensize'));
for i = 1:time
    if mod(i,3) == 0
        qx1surf = [];
        qz1surf = [];
        qx2surf = [];
        qz2surf = [];
        qx3surf = [];
        qz3surf = [];
        qx4surf = [];
        qz4surf = [];
        for l = 1:N
            qx1surf = [qx1surf,qx1(:,i)];            
            qz1surf = [qz1surf;linspace(0,0,50)];
            qx2surf = [qx2surf,qx2(:,i)-0.005];            
            qz2surf = [qz2surf;linspace(0,0,50)];
            qx3surf = [qx3surf,qx3(:,i)-0.01];            
            qz3surf = [qz3surf;linspace(0,0,50)];
            qx4surf = [qx4surf,qx4(:,i)-0.015];            
            qz4surf = [qz4surf;linspace(0,0,50)];
        end
        [mg11,mg12] = meshgrid(linspace(0,0.2,50),qy1(:,i));
        [mg21,mg22] = meshgrid(linspace(0,0.2,50),qy2(:,i));
        [mg31,mg32] = meshgrid(linspace(0,0.2,50),qy3(:,i));
        [mg41,mg42] = meshgrid(linspace(0,0.2,50),qy4(:,i));
        [mg51,mg52] = meshgrid(linspace(-0.35,0.65,500),linspace(-0.35,0.65,500));
        m1 = surf(mg11,mg12,qx1surf);
        hold on
        m2 = surf(mg21,mg22,qx2surf);
        m3 = surf(mg31,mg32,qx3surf);
        m4 = surf(mg41,mg42,qx4surf);
        m5 = surf(mg51,-0.301*ones(500,500),mg52);
        axis([-0.31,0.4,-0.31,0.4,-0.1,0.6])
        %axis([qz(1,i)-0.05,qz(1,i)+0.25,qy(1,i)-0.05,qy(1,i)+0.25,qx(1,i)-0.15,qx(1,i)+0.15,])
        view([80,170])
        grid on
        camroll(270)
        h = findobj('Type', 'surface');
        set(h(1), 'CData', original5, 'FaceColor', 'texturemap')
        set(h(2), 'CData', original1, 'FaceColor', 'texturemap')
        set(h(3), 'CData', original2, 'FaceColor', 'texturemap')
        set(h(4), 'CData', original3, 'FaceColor', 'texturemap')
        set(h(5), 'CData', original4, 'FaceColor', 'texturemap')
        set(m1, 'edgecolor','none')
        set(m2, 'edgecolor','none')
        set(m3, 'edgecolor','none')
        set(m4, 'edgecolor','none')
        set(m5, 'edgecolor','none')
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame);
        hold off
    end  
end

%%

