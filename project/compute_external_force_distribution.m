function e_dis = compute_external_force_distribution(jlist,N)
% the force will be exerted to the card at the first and the last node in the first 0.02 seconds
% after 0.02s the force will be canceled immediately.
time= size(jlist,2);
e_dis = zeros(3*N,time);
f = 1e-3;
for i = jlist
    if i< round(time/50)
       e_dis(1:3*3,i) = [f,0,-f,f,0,-f,f,0,-f];
       e_dis(3*(N-2)-2:3*N,i) = [f,0,f,f,0,-f,f,0,-f];
    end
end


