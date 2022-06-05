function [q_con,stack,thumb,n,stack_bound] = compute_stack(q_new)
% x coord is the function of y coord
% n is the power
n = 2;
global N
q_con = zeros(N,2);
for i = 1:N
    % define the constraint in discrete way
    q_con(i,:) = q_new(2*i-1:2*i)'+[-0.0003 0];
end
% use a polinomial equation to fit the constraint
% stack is a 3x1 vector, which represents 3 coefficient of a quatratic
% function
stack = polyfit(q_con(:,2),q_con(:,1),n);
thumb = q_con(1,2);
stack_bound = [q_con(1,2),q_con(N,2)];
end