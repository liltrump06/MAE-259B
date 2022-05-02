node0 = [0.1 0.1 0.1];
node1 = [0.2 0.2 0.1];
node2 = [0.2 0.1 0.2];
m1e = [0 0 1];
m2e = [0 0.7071 0.7071];
m1f = [0 0 1];
m2f = [0 0.7071 0.7071];
kappaBar = [0.1 0.1];
l_k = 0.1;
EI = 1e7;
[dF, dJ] = gradEb_hessEb(node0, node1, node2, ...
    m1e, m2e, m1f, m2f, kappaBar, l_k, EI);
theta_e = 0.5;
theta_f = 0.5;
refTwist = 0.2;
GJ = 3e6;
[F,J] = gradEt_hessEt(node0, node1, node2,theta_e, theta_f, refTwist,l_k, GJ)