function [Ff,Jf] = compute_Ff_Jf(Forceall,conmap,J,M_modinv,qd_old)
global N miu dt M
Ff = zeros(2*N,1);
Jf = zeros(2*N,2*N);
for i = 1:N
    frec = Forceall(2*i);
    v = (qd_old(2*i-1));
    if (conmap(i,:) == 1) && (abs(v)>1e-8) && sign(frec)>0
        sig = sign(v);
        Ff(2*i-1:2*i) = [-sig*frec*miu;0];
        Jfm = [-sig*miu*J(2*i,2*i-1),-sig*miu*(M(2*i,2*i)/(dt^2)+J(2*i,2*i));0,0];
        Jf(2*i-1:2*i,2*i-1:2*i) = dt*M_modinv(2*i-1:2*i,2*i-1:2*i)*Jfm;
    end
end
end