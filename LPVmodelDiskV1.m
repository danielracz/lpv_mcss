%Main script to run the example
%Requires the LPV core package  https://gitlab.com/tothrola/LPVcore
Ts=0.01;
Km=15.3145;
J=2.2*10^(-4);
m=0.1;
l=0.42/1000;
g=9.8;
tau=0.1; 


A_0=[-1/tau,0
      1     -1/tau]
A_1=[0,-m*g*l/J; 0,0]
B_0=[Km/tau;0]
B_1=[0;0]


%

np=1;
C_0=[0,1];
C_1=[0,0]
gamma=1.4*1/tau;
Gamma=gamma-5;
disp("Is max(norm(A_1,2)^2, (n_p+1)max(norm(C_1,2)^2,norm(C_2,2)^2) < Gamma")
disp(max([norm(A_1,2)^2, (np+1)*max([norm(C_0,2)^2,norm(C_1,2)^2])]) < Gamma)
disp("Is exp(A_0t) < e^{\gamma/2 t}")
disp(norm(expm(A_0),2) < exp(-gamma/2))
lambda=-(Gamma-gamma)*0.24

K_C=1
K_B=norm(B_0,2);


H2norm0=(np+1)*Gamma*K_B^2*(1/(gamma-lambda+Gamma)) %Estimate of H_2 norm using corrected version of Lemma 5.3


Q=eye(2)

P1=A_0'*Q+Q*A_0+A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q
P2=A_0'*Q+Q*A_0+0*A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q

H_2norm2=K_B^2*(np+1)*norm(Q)

[H_2norm22,Q22]=checkConditions(A_0,A_1,B_0,B_1,C_0,C_1,lambda)

np=2;
         
prevstream = RandStream.setGlobalStream(RandStream('twister','Seed',0)); %stream


th = preal('p','dt','Range', [-1/(2*pi-pi/2), 1]) % Define DT scheduling variable
A=eye(2)+Ts*[-1/tau, -(m*g*l/J)*th; 1 -1/tau];
B=Ts*[Km/tau; 0];
C=[0,1];
D=0;


Add(:,:,1)=A.peval(0);
Add(:,:,2)=A.peval(1)-Add(:,:,1);
Bdd(:,:,1)=B;
Bdd(:,:,2)=zeros(2,1);
Cdd(:,:,1)=C;
Cdd(:,:,2)=C;





[sysd_t,At_0,At_1,Bt_0,Bt_1,Ct_0,Ct_1]=converToCT(Add(:,:,1),Add(:,:,2),Bdd(:,:,1),Bdd(:,:,2),Cdd(:,:,1),Cdd(:,:,2),Ts);


[H2inf_true2,Qtrue2]=checkConditions(At_0,At_1,Bt_0,Bt_1,Ct_0,Ct_1,lambda);


sys_ss=LPVcore.lpvss(A,B,C,D,Ts);

pp = preal('p','ct','Range', [-1/(2*pi-pi/2), 1]) % Define DT scheduling variable
sys_ss_cont=LPVcore.lpvss(A_0+A_1*pp,B_0,C_0,0);

N=45;
N0=30;
%M=1000;


%NumberOfTs=[100,1000,5*10^3,10^4,2*10^4,5*10^4];
NumberOfTs=100:50:10*10^3; %,10^4,2*10^4,5*10^4];
%NumberOfTs=[100,200,300,500,550,600,650,700,2000,3000,4000,5000,10^4];
%NumberOfTs=[100,200,100,2000];
H2matrix=[] %H_2norm22];
ge_errs=[];
emp_errs=[];
ge_errs_1=[];
emp_errs_1=[];
maxoutputs=[];
M=max(NumberOfTs);

ParamEstError=[];
Bounds=[];
delta=0.05;

[ytrain,ytrain_cont,utrain,ptrain,Ybound,Ubound,CTimes] = SimulateTrueSystem(sys_ss,sys_ss_cont,2*M,N,Ts);sig(1,1)=1;
psig(1,1)=1;
psig(2,1)=var(ptrain(2,end,:));     
Ddd=0;


OutputNoises=0.05*randn(size(ytrain,1),size(ytrain,2),size(ytrain,3)); 

for j=1:size(NumberOfTs,2)



n=2;           % model order
p_window=3;   


nbpoints=NumberOfTs(j);

ytrain_local_clean=ytrain(:,N0+1:end,1:nbpoints);
ytrain_local=ytrain_local_clean+OutputNoises(:,N0+1:end,1:nbpoints);
utrain_local=utrain(:,N0+1:end,1:nbpoints);
ptrain_local=ptrain(:,N0+1:end,1:nbpoints); 


fprintf("Number of data points, %d",nbpoints);


[A0est,A1est,B0est]=lpvARXes_TS2(ytrain_local,utrain_local,ptrain_local,Ts);

B1est=[0;0];
C0est=[0,1];
C1est=[0,0];
 

sysEstDt=LPVcore.lpvss(A0est*Ts+eye(2)+Ts*A1est*th,B0est*Ts,C0est,0,Ts);



[H2inf,Qest]=checkConditions(A0est,A1est,B0est,B1est,C0est,C1est,lambda);


H2matrix=[H2matrix,H2inf];

param_error=[]

%for i=1:np
    param_error=[param_error, norm(A_0-A0est,2)]; %max(norm(Asig_d_true(:,:,i),2),0.0001)]
    param_error=[param_error, norm(A_1-A1est,2)]; %max(norm(Bsig_d_true(:,:,i),2),0.0001)];
    param_error=[param_error, norm(B_0-B0est,2)]; %max(norm(Csig_d_true(:,:,i),2),0.0001)];
%end

ParamEstError=[ParamEstError, max(param_error)];

if (H2inf==-1)
  fprintf('Estimated model does not satisfy the condition, %f',M)
end



[ge_err,emp_err,maxsimoutput, ge_err_1,emp_err_1]=GenerateError(sysEstDt,nbpoints,utrain,ptrain,ytrain+OutputNoises,CTimes);

ge_errs=[ge_errs,ge_err];
emp_errs=[emp_errs,emp_err];
ge_errs_1=[ge_errs_1,ge_err_1];
emp_errs_1=[emp_errs_1,emp_err_1];
maxoutputs=[maxoutputs,maxsimoutput];

Bounds=[Bounds, 1/sqrt(nbpoints)*(2+sqrt(2*log(4/delta)))]
end

fprintf ("H2 norms");
H2matrix
fprintf("Parameter estimation error")
ParamEstError


c1=max(H2matrix);
Lu=Ubound;
c2=max(max(max(abs(ytrain))), max(maxoutputs));
Kl=2*c2;
c=2*Kl*max([Lu*(np)*c1,c2]);
Bounds_2=c*Bounds; 
Bounds_1=2*max([Lu*(np)*c1,max(max(max(abs(ytrain))))])*Bounds;

save("Data.mat","NumberOfTs", "ge_errs","Bounds", "Bounds_1", "Bounds_2", "ge_errs_1","emp_errs","emp_errs_1");

save("DataSecond.mat", "H2matrix","ParamEstError","utrain","ytrain","ptrain","OutputNoises","maxoutputs","ytrain_cont","utrain_local","ptrain_local");
PlotData("Data.mat");



return;


