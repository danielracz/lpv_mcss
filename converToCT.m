function [sysd,A0,A1,B0,B1,C0,C1]=converToCT(Ad0,Ad1,Bd0,Bd1,Cd0,Cd1,Ts)

A0=(Ad0-eye(size(Ad0,1)))/Ts;
A1=Ad1/Ts;
B0=Bd0/Ts;
B1=Bd1/Ts;
C0=Cd0;
C1=Cd1;
p  = preal('p','ct','Range', [-1/(2*pi-pi/2), 1]);

A=A0+A1*p;
B=B0+B1*p;
C=C0+C1*p;
D=0;

sysd=LPVcore.lpvss(A,B,C,D);

end