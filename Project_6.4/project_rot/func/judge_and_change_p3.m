function [conmap,recal] = judge_and_change_p3(velop, ...
    q_new,q_old,conmap,F)
global N dt ground
%% add constraint
recal = 0;
variance = 0.00001;
for i = 1:N
    ypos = q_new(2*i);
    if conmap(i) == 0
        if ypos < ground
            fprintf('constraint added at node %d caused by ground \n',i)
            conmap(i) =1;
            recal = 1;
        end
    end
    % calculate the unit normal vector at contact point
end

%% remove constraint
for i = 1:N
    ypos = q_new(2*i);
    if conmap(i) == 1
        frec = F(2*i-1:2*i);
        fnormal = dot(frec,[0,1]);
        if (fnormal <=0) && (ypos >= ground - 2e-6)
            fprintf('constraint deleted at node %d \n',i)
            conmap(i) = 0;
            recal = 1;
        end
    end
end

end


    
