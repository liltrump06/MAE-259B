function contact_standard(q_old,q_new,stack,thumb)
global N
%contact_ind = zeros(N,1);
n = size(stack,2)-1;
for k = 1:N
    q = q_new(4*k-3:4*k-1);
    for i = 0:n
        line = line + stack(n-i+1)*q(3).^i;
    end
    if line > q_new(1)
        
    end
end



