function [ytrain,ytrain_cont,utrain,ptrain,Ybound,Ubound,CTinterval] = SimulateTrueSystem(sys_ss,sys_ss_cont,M,N,Ts)
ytrain=zeros(1,N,M);
%ytrain_real=zeros(1,N,M);
utrain=zeros(1,N,M);
ptrain=ones(2,N,M);

Ybounds=zeros(1,M);
Ubounds=zeros(1,M);
CTinterval=0:Ts:Ts*(N-1);
ytrain_cont=zeros(1,N,M);
for i=1:M
    %ut=idinput(N,'rgs');
    %p_mag=0.25;
    %pt=(1-p_mag)+p_mag*idinput(N,'rbs',[0,0.05]);
    pt=ones(2,N);
    pt(2,:)=3*rand(1,N)-1.5;
    ut=30*(rand(N,1)-0.5);
    ysim=lsim(sys_ss,pt(2,:)',ut); 
    %ysim_real=lsim(sys_real_true,pt(2,:)',ut);
    ysim_cont=ysim; %lsim(sys_ss_cont,pt(2,:)',ut,CTinterval);
    ytrain(:,:,i)=ysim';
    %ytrain_real(:,:,i)=ysim_real;
    ytrain_cont(:,:,i)=ysim_cont;
    utrain(:,:,i)=ut';
    ptrain(:,:,i)=pt;
    Ybounds(1,i)=norm(ysim_cont,Inf);
    Ubounds(1,i)=norm(ut,2)*Ts;
end
Ybound=max(Ybounds);
Ubound=max(Ubounds);


end