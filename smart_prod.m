function [res]=smart_prod(lambda)
   slambda=sort(lambda);

   N=size(lambda,2)
   i=floor(N/2);
  
   res=0;
   i
   for k=1:i
       %resold=res;
       res=res+log(slambda(1,k))+log(slambda(1,N-k+1));
       disp('k')
       k
       disp('N-k+1')
       N-k+1
       if abs(res) < 1e-200
            disp('k')
            k
            disp('first')
            slambda(1,k)
            disp('second')
            slambda(1,N-k+1)
            disp('resold')
            resold
            res=exp(res)
            return
        end    
        if abs(res) > 1e300
            disp('k')
            k
            disp('first')
            slambda(1,k)
            disp('second')
            slambda(1,N-k+1)
            disp('resold')
            resold
            res=exp(res)
            return
        end
   end
   disp('End')
   res
   res=exp(res)
  end