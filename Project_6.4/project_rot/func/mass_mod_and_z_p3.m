function [imposeacc,M_modinv,nlist,interlist] = mass_mod_and_z_p3(q_new,M,conmap,velop)
global N ground dt
imposeacc = zeros(2*N,1);
M_modinv  = zeros(2*N,1);
nlist = zeros(N,2);
interlist = zeros(N,2);
for i =1:N
    if conmap(i) == 1
        xpos = q_new(2*i-1);
        ypos = q_new(2*i);
        vpoint = velop(2*i-1:2*i);
        y0 = ground;
        x0 = xpos;
        normvec = [0;1];
        if normvec(1) < 0
            normvec = -normvec;
        end

        nlist(i,:) = normvec;
        imposeacc(2*i-1:2*i) = ([x0;y0] - [xpos;ypos]) / dt;
        %imposeacc(2*i-1:2*i) = ([x0;y0] - [xpos;ypos]) / dt - normvec * dot(vpoint, normvec);
        S = Sformod(normvec);
        M_modinv(2*i-1:2*i,2*i-1:2*i) = M(2*i-1:2*i,2*i-1:2*i)\S;
        interlist(i,1) = x0;
        interlist(i,2) = y0;
    else
        imposeacc(2*i-1:2*i) = 0;
        M_modinv(2*i-1:2*i,2*i-1:2*i) = inv(M(2*i-1:2*i,2*i-1:2*i));
    end

end

end


