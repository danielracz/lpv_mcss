%Markov chain
%Define the transition probability matrix
Ntot=5e4
np=2
prob.s11=1/4; prob.s12=1-prob.s11; prob.s21=1/3; prob.s22=1-prob.s21;
%prob.s11=1/3; prob.s12=1-prob.s11; prob.s21=2/3; prob.s22=1-prob.s21
%prob.s11=1/2; prob.s12=1-prob.s11; prob.s21=1/2; prob.s22=1-prob.s21
P = [prob.s11, prob.s12;  % Probability of transitioning from state 1 to state 1 and state 2
     prob.s21, prob.s22]; % Probability of transitioning from state 2 to state 1 and state 2
% CrXeate a Discrete-Time Markov Chain object
mc = dtmc(P);
% Generate a sequence of states
numSteps = Ntot;
qq = (simulate(mc, numSteps))';

p=zeros(np,Ntot);
for t=1:Ntot+1
    for j=1:np
        if qq(1,t)==j
            p(j,t)=1;
        else
            p(j,t)=0;
        end
    end
end


%p=zeros(np,Ntot)

psig = zeros(np,1);
% psig(1,:)=0.5;
for i=1:np
    psig(i,1) = 1/Ntot*sum((p(i,:).*p(i,:)));     
    %var(p(i,:).*sqrt(lambda)) + (mean(p(i,:).*lambda))^2;
end
% psig(2,1)=0.5;
% lambda=ones(1,Ntot);

lambda=ones(1,Ntot);
for t=2:Ntot
    for i=1:np
        if(qq(1,t)==i)
            for j=1:np
                if(qq(1,t-1)==j)
                    pij=1/prob.(['s',sprintf('%.0f',j),sprintf('%.0f',i)]);
                    lambda(t)=(1/prob.(['s',sprintf('%.0f',j),sprintf('%.0f',i)]));

                end
            end
        end
    end
end

lambda_train=lambda(floor(Ntot/2):Ntot)*;
qqtrain=qq(floor(Ntot/2):Ntot);

ptrain=p(floor(Ntot/2):Ntot);

% zz=lambda_train(1:50)
% zz(zz >1)
% prod(zz(zz >1))
% prod(zz(zz < 1))


z1=(qqtrain==2);
z1=z1(2:end-1);
z2=(qqtrain==1);
z3=z2(1:end-2);
z2=z2(3:end);

mean(z1.*z2.*z3.*lambda_train(3:end).*lambda_train(2:end-1).*lambda_train(1:end-2))/mean(z3.*lambda_train(1:end-2).*lambda_train(2:end-1))


mean(z3.*lambda_train(3:end).*lambda_train(2:end-1).*lambda_train(1:end-2))

T0=50
z4=(qqtrain==2);
zlambda=ones(1,size(lambda_train,2)-T0);
for i=1:size(lambda_train,2)-T0
     zlambda(1,i)=exp(sum(log(lambda_train(i:T0+i))));
end    
disp('Mean past')
mean(z4(T0+1:end).*zlambda)


T0=500
z4=(qqtrain==2);
zlambda=ones(1,size(lambda_train,2)-T0);
for i=1:size(lambda_train,2)-T0
     zlambda(1,i)=exp(sum(log(lambda_train(i:T0+i))));
end    
disp('Mean')
mean(z4(1:end-T0).*zlambda)

disp("Product")
zz=lambda_train(5e3:9e3);
zz(zz >1);
prod(zz(zz >1))
prod(zz(zz < 1)) 
