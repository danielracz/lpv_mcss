function [lpvsys,A_0,A_1,B_0,B_1,C_0,C_1]=CreateLPV(J,Km,g,m,tau,l)
A_0=[-1/tau,0
      1     -1/tau]
A_1=[0,-m*g*l/J; 0,0]
B_0=[Km/tau;0]
B_1=[0;0]
%
np=1;
C_0=[0,1];
C_1=[0,0];
p  = preal('p','ct','Range', [-1/(2*pi-pi/2), 1]);
A=A_0+A_1*p;
B=B_0+p*B_1;
C=C_0+p*C_1;
lpvsys=LPVcore.lpvss(A,B,C,0);
end