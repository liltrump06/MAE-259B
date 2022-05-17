function [q_con,stack,thumb,n] = compute_stack(q_new)
% x coord is the function of z coord
% n is the power
n = 2;
global N
q_con = zeros(N,2);
for i = 1:N
    q_con(:,i) = q_new(2*i-1:2*i)+[0.0003 0]';
end
stack = polyfit(q_con(:,2),q_con(:,1),n);
thumb = q_con(1,2);
end