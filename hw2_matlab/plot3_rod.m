clc; clear; close all
load('hw2p2_data.mat')
N = (size(q_whole,1)+1)/4;
time = size(q_whole,2);
qplot = zeros(N,3,time);
for i = 1:N
    for t = 1:time
        qplot(i,:,t) = q_whole(4*i-3:4*i-1,t);
    end
end
v = VideoWriter('Output2.avi');
open(v);
figure(1)
set(gcf,'outerposition',get(0,'screensize'));
for t = 1:size(q_whole,2)
    plot3(qplot(:,1,t),qplot(:,2,t),qplot(:,3,t),'-o')
    axis([-0.05,0.05,-0.05,0.05,-0.1,0.05])
    view([10,0])

    frame = getframe(gcf);
    writeVideo(v,frame);
    drawnow
    
end

