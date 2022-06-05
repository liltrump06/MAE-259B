function Ff = compute_Ff(Forceall,conmap,qd_old)
global N miu 
Ff = zeros(2*N,1);
Jf = zeros(2*N,2*N);
for i = 1:N
    frec = Forceall(2*i);
    v = (qd_old(2*i-1));
    if (conmap(i,:) == 1) && (abs(v)>1e-8) && sign(frec)>0
        sig = sign(v);
        Ff(2*i-1:2*i) = [-sig*frec*miu;0];
    end
end