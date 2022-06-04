global N l r0 rou E I A rou_al dampratio
dampratio = 1e-7;
N = 50;
l = 0.14605; % length
dl = l/(N-1);
r0 = 1.9074e-03/2; % radius of the rod
rou = 936; % density kg/m3
E = 1e8; % Young's Modulus
I = pi/4*(r0^4);
A = pi*r0^2;
Rc = [cos(-pi/6),-sin(-pi/6);sin(-pi/6),cos(-pi/6)];
m  = pi*(r0^2)*l*rou/(N-1);
q = zeros(2*N,1); 
miu = 0.05;
