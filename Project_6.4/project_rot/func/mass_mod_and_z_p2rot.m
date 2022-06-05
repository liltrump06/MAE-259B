function [z,M_modinv,nlist,interlist] = mass_mod_and_z_p2rot(q_new,M,conmap,velop)
global N stack dt thumb Rc
z = zeros(2*N,1);
M_modinv  = zeros(2*N,1);
nlist = zeros(N,2);
interlist = zeros(N,2);
qr_new = q_new;
% convey the real coordinates to the original frame
for i = 1:N
    q_new(2*i-1:2*i) = Rc' * qr_new(2*i-1:2*i);
end
for i =1:N
    if conmap(i,1) == 1
        if conmap(i,2) == 0
            xpos = q_new(2*i-1);
            ypos = q_new(2*i);
            vpoint = velop(2*i-1:2*i);
            % 3 coefficients
            a = stack(1);
            b = stack(2);
            c = stack(3);
            % 3 possible y positions for the q_correct point.
            rlist = roots([2*a^2,3*a*b,2*a*c+b^2-2*a*xpos+1,-b*xpos+c*b-ypos]);
            length = 100;
            % sift for the only possible y coordinate
            for l = 1:size(rlist,1)
                m = rlist(l);
                if imag(m)==0
                    if abs(m-ypos) < length
                        length = abs(m-ypos);
                        y0 = m;
                    end
                end
            end
            % compute the x position of q_correct
            x0 = y0^2*stack(1)+y0*stack(2)+stack(3);
            % compute the unit constraint direction vector in real world frame
            normvec = Rc*[x0-xpos;y0-ypos]/norm(Rc*[x0-xpos;y0-ypos]);
            if normvec(1) < 0
                normvec = -normvec;
            end
            nlist(i,:) = normvec;
            % compute the compensated velocity vector
            z(2*i-1:2*i) = Rc*([x0;y0] - [xpos;ypos]) / dt;
            % compute modified inverted mass matrix
            S = Sformod(normvec);
            M_modinv(2*i-1:2*i,2*i-1:2*i) = M(2*i-1:2*i,2*i-1:2*i)\S;
            interlist(i,1) = x0;
            interlist(i,2) = y0;
        elseif conmap(i,2) == 1
            xpos = q_new(2*i-1);
            ypos = q_new(2*i);
            vpoint = velop(2*i-1:2*i);
            % unit constraint direction vector is a constant vector for
            % thumb constraint
            normvec =Rc*[0,1]';
            nlist(i,:) = normvec;
            x0 = xpos;
            y0 = thumb;
            z(2*i-1:2*i) = Rc*([x0;y0] - [xpos;ypos]) / dt - normvec * dot(vpoint, normvec);
            S = Sformod(normvec);
            M_modinv(2*i-1:2*i,2*i-1:2*i) = M(2*i-1:2*i,2*i-1:2*i)\S;
        end
    else 
        z(2*i-1:2*i) = 0;
        M_modinv(2*i-1:2*i,2*i-1:2*i) = inv(M(2*i-1:2*i,2*i-1:2*i));

    end

end


