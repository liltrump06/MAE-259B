N = 5;
q_new = [ 2e-02,0,0,0,6.306e-03,1.898e-02,0,0,-1.602e-02,1.197e-02,...
    0,0,-1.641e-02,-1.143e-02,0,0,5.673e-03,-1.918e-02,0;2e-02,0,0,0,6.306e-03,1.898e-02,0,-2.006e-01,-1.463e-02,1.281e-02,-8.462e-03,2.191e-01,...
    -1.441e-02,-8.443e-03,-1.827e-02,4.726e-01,7.321e-03,-1.712e-02,-1.796e-02;2e-02,0,0,0,6.306e-03,1.898e-02,0,-3.908e-01,-1.119e-02,1.474e-02,-1.497e-02,3.572e-01,...
    -9.813e-03,3.396e-04,-3.337e-02,9.616e-01,1.117e-02,-9.421e-03,-3.688e-02]';
u11 = [-8.110e-01,-5.851e-01,0];
u1 = zeros(N-1,3);
u2 = zeros(N-1,3);
tk_old = zeros(N-1,3,3);
tk_new = zeros(N-1,3,3);
a1_old = zeros(N-1,3,3);
a1_new = zeros(N-1,3,3);
a2_old = zeros(N-1,3,3);
a2_new = zeros(N-1,3,3);
m1 = zeros(N-1,3);
m2 = zeros(N-1,3);
mkref = zeros(N-2,3);
for j = 1:3
    for m = 1:N-1
        tk1 =(q_new(4*m+1:4*m+3,j) - q_new(4*m-3:4*m-1,j))/norm((q_new(4*m+1:4*m+3,j) - q_new(4*m-3:4*m-1,j)));
        if j == 1
            if m == 1
                u1(m,:,j) = u11;
            else
                u1(m,:,j) = parallel_transport(u1(m-1,:,j)',tk_new(m-1,:,j)',tk1);
            end
            a1_new(m,:,j) = u1(m,:,j); 
        else
            if m == 1
                u1(m,:,j)= u11;
            else
                u1(m,:,j)= parallel_transport(u1(m-1,:,j)',tk_new(m-1,:,j)',tk1);
            end
            a1_new(m,:,j)= parallel_transport(a1_new(m,:,j-1)',tk_new(m,:,j-1)',tk1);
        end
        a2_new(m,:,j)= cross(tk1,a1_new(m,:,j));
        u2(m,:,j)= cross(tk1,u1(m,:,j));
        tk_new(m,:,j)= tk1;
        m1(m,:,j)= cos(q_new(4*m,j))*a1_new(m,:,j)+sin(q_new(4*m,j))*a2_new(m,:,j);
        m2(m,:,j)= cos(q_new(4*m,j))*a2_new(m,:,j)-sin(q_new(4*m,j))*a1_new(m,:,j);
    end
    for k = 1:N-2
        mkref(k,j) = signedAngle(parallel_transport(a1_new(k,:,j)',tk_new(k,:,j)', ...
        tk_new(k+1,:,j)'),a1_new(k+1,:,j)',tk_new(k,:,j)');
    end
end
