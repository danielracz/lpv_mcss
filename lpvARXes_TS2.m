function [A_0,A_1,B_0]=lpvARXes_TS2(ytrain,utrain,ptrain,Ts)

th=preal('p','dt');
th1=pshift(th,-1);
th2=pshift(th,-2);
am1=1+th+th1+th2;

am2=am1;
bm0=am1;
bm1=am1;
bm2=am1;
template_sys = LPVcore.lpvidpoly({1, am1, am2}, {bm0, bm1, bm2});

% perform ARX identification with the specified model structure
options = lpvarxOptions('Display','on');


Y=reshape(ytrain(1,end,:),size(ytrain,3),1);
y_1=reshape(ytrain(1,end-1,:),size(ytrain,3),1);
y_2=reshape(ytrain(1,end-2,:),size(ytrain,3),1);
p_2=reshape(ptrain(2,end-2,:),size(ptrain,3),1);
u_2=reshape(utrain(1,end-2,:),size(utrain,3),1);
X=[y_1,y_2, y_2.*p_2, u_2];
%X=[y_1,y_2,u_2];
%Y=Y+-1.8709*y_2.*p_2;



theta=lscov(X,Y);

a1=(theta(1)/2-1)/Ts;
%a1=min(-4.9,max(-10,a1));
a2=theta(3)/Ts^2;
%a2=max(-2,min(-2,a2));
b=theta(3)/Ts^2;

A_0=[a1,0
     1,a1];

A_1=[0,a2
     0 0];

B_0=[b;0]

%y=reshape(ytrain, size(ytrain,1)*size(ytrain,2)*size(ytrain,3),1);

%u=reshape(utrain, size(utrain,1)*size(utrain,2)*size(utrain,3),1);
%p=reshape(ptrain(2,:,:), size(utrain,2)*size(utrain,3),1);

%data_est = lpviddata(y,p,u,Ts);  


%m_arx = lpvarx(data_est, template_sys, options);

%alpha1=m_arx.A.Mat{2}.matrices(1);

end
