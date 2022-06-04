function [conmap,velop,recal] = judge_and_change_p2rot(velop, ...
    q_new,q_old,conmap,F)
global N dt stack stack_bound thumb Rc
%% add constraint
for i = 1:N
    q_new(2*i-1:2*i) = Rc' * q_new(2*i-1:2*i);
end
recal = 0;
variance = 0.00001;
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    line = stack(1)*ypos^2+stack(2)*ypos+stack(3);
    if conmap(i,1) == 0
        if ypos < thumb && xpos < 0.02
            fprintf('constraint added at node %d caused by thumb \n',i)
            conmap(i,1) =1;
            conmap(i,2) =1;
            velop(2*i-1,:) = (q_new(2*i-1)-q_old(2*i-1))/dt;
            velop(2*i,:) = (q_new(2*i)-q_old(2*i))/dt;
            recal = 1;
        elseif (ypos > stack_bound(1)-variance && ypos < stack_bound(2)+variance && xpos < line)
            fprintf('constraint added at node %d, caused by stack\n',i)
            conmap(i,1) =1;
            conmap(i,2) =0;
            velop(2*i-1,:) = (q_new(2*i-1)-q_old(2*i-1))/dt;
            velop(2*i,:) = (q_new(2*i)-q_old(2*i))/dt;
            recal = 1;
        end
    end
    % calculate the unit normal vector at contact point
end

%% remove constraint
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    line = stack(1)*ypos^2+stack(2)*ypos+stack(3);
    if conmap(i,1) == 1
        frec = F(2*i-1:2*i);
        if conmap(i,2) == 1 % thumb cons 
            fnormal = dot(frec,[0,1]);
            if (fnormal <=0) && (ypos >= thumb - 2e-6 || xpos > 0.02)
                fprintf('constraint deleted at node %d \n',i)
                conmap(i,1) = 0;
                recal = 1;
            end
        elseif conmap(i,2) == 0 % stack cons
            fnormal = dot(frec,nstack(ypos,xpos));
            if (fnormal <=0) && (xpos >= line - 2e-6 || (ypos < stack_bound(1)-variance || ypos > stack_bound(2)+variance))
                fprintf('constraint deleted at node %d \n',i)
                conmap(i,1) = 0;
                recal = 1;
            end
        end
    end

end


    
