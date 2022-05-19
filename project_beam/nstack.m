function normvec = nstack(ypos,xpos)
global stack
rlist = roots([2*stack(1)^2,3*stack(1)*stack(2),2*stack(1)*stack(3)-stack(2)^2+2*stack(1)*xpos-1,stack(2)*xpos-stack(3)*stack(2)+ypos]);
for i = 1:size(rlist,1)
    m = rlist(i);
    length = 100;
    if imag(m)==0
        if abs(m-ypos) < length
            length = abs(m-ypos);
            y0 = m;
        end
    end
end      
x0 = y0^2*stack(1)+y0*stack(2)+stack(3);
normvec = [x0-xpos;y0-ypos]/norm([x0-xpos;y0-ypos]);
