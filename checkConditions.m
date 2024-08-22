function [H2inf,Q]=checkConditions(A0,A1,B0,B1,C0,C1,lambda)

%P1=A_0'*Q+Q*A_0+A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q


%P2=A_0'*Q+Q*A_0+0*A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q

setlmis([]);
X = lmivar(1,[2 1]);

epsilon=0.01;
lmiterm([1 1 1 X],1,A0,'s');
lmiterm([1 1 1 X],A1',A1);
lmiterm([1 1 1 0],C0'*C0);
lmiterm([1 1 1 0],C1'*C1);
lmiterm([1,1,1,X],lambda,1);
lmiterm([1,1,1,0],epsilon*eye(2))
lmiterm([1,2,2,X],-1,1);

lmis = getlmis;

n = decnbr(lmis); 
c = zeros(n,1);

for j=1:n 
	Xj = defcx(lmis,j,X);
	c(j) = B0'*Xj*B0+B1'*Xj*B1;
end

options = [1e-5,0,0,0,0] ;
[copt,xopt] = mincx(lmis,c,options);

%[tmin,xfeas] = feasp(lmis);

H2inf=-1;
Q=zeros(2,2);
if size(xopt,1) ~= 0

    Q=dec2mat(lmis,xopt,X);
    disp('eiegnevalues Q')
    %eig(Q)
    %disp('Q')
    Q;
    H2inf=trace(B0'*Q*B0+B1'*Q*B1);
    %disp('Copt')
    %disp(copt);
    %P2=A0'*Q+Q*A0+A1'*Q*A1+C1'*C1+C0'*C0+lambda*Q
    %disp('eige P2')
    %eig(P2)
end
end
