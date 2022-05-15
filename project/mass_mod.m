function M_mod =mass_mod(M,plist,qlist)
%%
% M [4N-1,4N-1]
% plist [3,N]
% qlist [3,N]
global N
for i = 1:N
    M1 = M(4*i-3:4*i-1,4*i-3:4*i-1);
    p1 = plist(:,i);
    q1 = qlist(:,i);
    Sinv = 1\(eye(3)-p1*p1'-q1*q1');
    M(4*i-3:4*i-1,4*i-3:4*i-1) = M1*Sinv;
end
M_mod = M; 
    