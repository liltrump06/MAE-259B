function [qx,qy,qz] = make_card_surf(qall,N)
for i = 1:N
    qx(i,:) = qall(2*i-1,:);
    qy(i,:) = qall(2*i,:);
end
qz = zeros(size(qy));