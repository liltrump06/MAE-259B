function [F,Jaco]= computeF_J(N,q,l_k,EA,m1,m2,GJ,EI,kappalist,mkref,g,M)
    F = zeros(4*N-1,1);
    Jaco = zeros(4*N-1,4*N-1);
    for k = 1:N-1
        ind = [4*k-3:4*k-1,4*k+1:4*k+3];
        [dF,dJ] = gradEs_hessEs(q(4*k-3:4*k-1)',q(4*k+1:4*k+3)',l_k,EA);
        F(ind)= F(ind)+dF;
        Jaco(ind,ind) = Jaco(ind,ind)+dJ;
    end
    for k = 1:N-2
        ind = 4*k-3:4*k+7;
        node0 = q(4*k-3:4*k-1)';
        node1 = q(4*k+1:4*k+3)';
        node2 = q(4*k+5:4*k+7)';
        [dF,dJ] = gradEb_hessEb(node0,node1,node2,m1(k,:),m2(k,:),m1(k+1,:),m2(k+1,:),kappalist(k,:),l_k,EI);
        F(ind) = F(ind)+dF;
        Jaco(ind,ind) = Jaco(ind,ind)+dJ;
        
        [dF,dJ] = gradEt_hessEt(node0,node1, ...
            node2,q(4*k),q(4*k+4),mkref(k),l_k,GJ);        
        F(ind)= F(ind)+dF;
        Jaco(ind,ind) = Jaco(ind,ind)+dJ;
    end
    for k = 1:N
        F(4*k-3:4*k-1) = F(4*k-3:4*k-1)- M(4*k-3:4*k-1,4*k-3:4*k-1)*g';
    end
    
end


