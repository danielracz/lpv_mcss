function [gen_error, emp_error, maxsimoutput, gen_error_1,emp_error_1]=GenerateError(sys_est, end_training,u,p,y,t)


gen_error=0;
emp_error=0;
gen_error_1=0;
emp_error_1=0;
maxsimoutput=0;
for j=1:size(u,3)
    ysim_est=lsim(sys_est,p(2,:,j)',u(:,:,j)'); %t);
    

    if j <= end_training
        emp_error=emp_error+norm(y(:,:,j)'-ysim_est,2)^2;
        emp_error_1=emp_error_1+norm(y(:,:,j)'-ysim_est,2);
    else
        gen_error=gen_error+norm(y(:,:,j)'-ysim_est,2)^2;
        gen_error_1=gen_error_1+norm(y(:,:,j)'-ysim_est,2);
    end

    if maxsimoutput < max(abs(ysim_est))
        maxsimoutput= max(abs(ysim_est));
    end

end

gen_error=gen_error/(size(u,3)-end_training);
emp_error=emp_error/end_training;
gen_error_1=gen_error_1/(size(u,3)-end_training);
emp_error_1=emp_error_1/end_training;
end