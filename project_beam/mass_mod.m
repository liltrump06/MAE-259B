function mass_mod_and_z(q_new,q_old,M,conmap,velop)
global N stack
for i =1:N
    k = (q_new(2*i-1)-q_old(2*i-1))/(q_new(2*i)-q_old(2*i));
    b = q_new(2*i-1);
    vp = velop(i,:);
    vp_normal = dot(vp,nstack(ypos))*nstack(ypos);


