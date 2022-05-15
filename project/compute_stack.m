
%[q_con,stack,thumb,n] = compute_stack1(q_new);
%km = linspace(0.0242,0.1461,50);
%line = 0;
%for i = 0:n
%    line = line + stack(n-i+1)*km.^i;
%end
%plot(q_con(3,:),q_con(1,:),'ro',km,line,'b-')


function [q_con,stack,thumb,n] = compute_stack(q_new)
% x coord is the function of z coord
% n is the power
n = 7;
global N
q_con = zeros(3,N);
for i = 1:N
    q_con(:,i) = q_new(4*i-3:4*i-1)+[0.0003 0 0]';
end
stack = polyfit(q_con(3,:),q_con(1,:),n);
thumb = q_con(3,1);
end

