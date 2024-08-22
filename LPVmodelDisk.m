Ts=0.1;
Km=15.3145*0.03;
J=2.2*10^(-4);
m=0.07;
l=0.42/1000;
g=9.8;
tau=0.5971*0.3;

%Continuous-time
%x(t)=[dtheta(t);theta(t)]
%p(t)=sint(theta(t))/theta(t))
%ddot theta(t)=-1/tau theta-m*g*l/Jptheta-Km/tau u, u=-1/Km theta+v
A_0=[-1/tau,0
      1     -1/tau]
A_1=[0,-m*g*l/J; 0,0]
B_0=[Km/tau;0]
B_1=[0;0]

%Disrete-time x(t)=[theta(kT),theta((k-1)T)]
Adt_0=[Ts/tau-2, 1-Ts/tau
     1       , 0];

A1dt=[0, m*g*l*Ts^2/J,
    0,0];

Bd0=[Km*Ts^2/tau;0];
Bd1=[0:0];

N = 5000;
u = sin((1:N)*Ts*2)';
p_mag=0.25;
p=(1-p_mag)+p_mag*cos((1:N)*Ts)';


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
lambda=-(Gamma-gamma)*0.3

%H2 norm 
K_C=1
K_B=norm(B_0,2);


H2norm=(np+1)*Gamma*K_B^2*(1/(gamma-lambda+Gamma)) %Estimate of H_2 norm using corrected version of Lemma 5.3


Q=eye(2)

P1=A_0'*Q+Q*A_0+A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q
P2=A_0'*Q+Q*A_0+0*A_1'*Q*A_1+C_1'*C_1+C_0'*C_0+lambda*Q

H_2norm2=K_B^2*(np+1)*norm(Q) %estimate of H_2 norm using Lemma 4.3

eig(P1)
eig(P2)


Ts = 0.1;              % Sampling time
p  = preal('p','dt','Range', [-1/(2*pi-pi/2), 1]);

p2= pmatrix(1, ...
    'SchedulingTimeMap',timemap([-2],'dt','Name',{'p'}));
                            % Dynamic dependence, shifted back 2 samples
 
% Alternatives:
p2= preal('p','dt','Dynamic',-2);  % Direct definition via preal
p2= pshift(p,-2);                     % By the Time Shift operator
                            
                       
a0=1;                       % Coefficients of the output side         
a1=Ts/tau-2;               
a2=1-Ts/tau+(m*g*l)/J*Ts^2*p2;  % Coefficient with dynamic dependence;
 
b2=Km/tau*Ts^2;   

th = preal('p','dt','Range', [-1/(2*pi-pi/2), 1]) % Define DT scheduling variable
A=eye(2)+Ts*[-1/tau, -(m*g*l/J)*th; 1 0];
B=Ts*[Km/tau; 0];
C=[0,1];
D=0;

Ts=0.01
N = 500;
M=1000;
sys_ss=LPVcore.lpvss(A,B,C,D,Ts);
ytrain=[];
utrain=[];
ptrain=[]
for i=1:M
    u = sin((1:N)*Ts*2)';
    p_mag=0.25;
    p=(1-p_mag)+p_mag*cos((1:N)*Ts)'; 
    ysim=lsim(sys_ss,p,u); 
    ytrain=[ytrain;ysim];
    utrain=[utrain;u];
    ptrain=[ptrain;p];    
end

DATA=lpviddata(ytrain,ptrain,utrain,Ts);

alpha=0.8;
Ainit=eye(2)+Ts*[-1/tau, -alpha*(m*g*l/J)*th; 1 0];
sys_init=LPVcore.lpvidss(Ainit,B,C,D, 'innovation', zeros(2,1), [],1,Ts);

SYS_Learned=lpvssest(DATA,sys_init,lpvssestOptions('Display','on','InitialState','zero'));

compare(DATA,SYS_Learned, struct('InitialCondition', 0));