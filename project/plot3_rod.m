clc; clear; close all
load('project.mat')
N = (size(q_whole,1)+1)/4;
time = size(q_whole,2);
qplot = zeros(N,3,time);
for i = 1:N
    for t = 1:time
        qplot(i,:,t) = q_whole(4*i-3:4*i-1,t);
    end
end
qcard = qplot;
q_para = qplot;
for i = 1:5
    q_para(:,2,:) = qplot(:,2,:)+0.01*i;
    qcard = cat(1,qcard,q_para);
end
%%
d = 0:0.05:2;
v = VideoWriter('Output2.avi');
open(v);
figure(2)
set(gcf,'outerposition',get(0,'screensize'));
for t = 1:size(q_whole,2)
    if mod(t,10) ==0
        scatter3(qcard(:,1,t),qcard(:,2,t),qcard(:,3,t),'o')
        axis([-1,0.5,-0.5,0.5,-0.5,0.5])
        view([10,10])
        frame = getframe(gcf);
        writeVideo(v,frame);
        drawnow
    end
end

