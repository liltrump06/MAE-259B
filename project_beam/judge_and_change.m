function [conmap,velop,recal] = judge_and_change(velop, ...
    q_new,q_old,conmap,F)
global N dt stack polyn
%% add constraint
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    for m = 0:polyn
        line = line + stack(n-i+1)*ypos.^i;
    end
    if conmap(i) == 0 && xpos < line
        fprintf('constraint added at node %d',i)
        conmap(i) =1;
        velop(i,:) = (q_new(2*i-1:2*i)-q_old(2*i-1:2*i))/dt;

        recal = 1;
        ny = -(2*stack(1)*ypos+stack(2));
        nx = 1;
        nlist(i,:) = [nx;ny]/norm([nx;ny]);
    end
    % calculate the unit normal vector at contact point
end

%% remove constraint
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    for m = 0:polyn
        line = line + stack(n-i+1)*ypos.^i;
    end
    if conmap(i) == 1
        frec = F(2*i-1:2*i);
        fnormal = dot(frec,nstack(ypos,stack));
        if fnormal <=0 && xpos >= line - 2e-6
            fprintf('constraint deleted at node %d',i)
            conmap(i) = 0;
            recal = 1;
        end
    end
end


    
