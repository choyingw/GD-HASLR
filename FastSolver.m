function [X,L] = FastSolver(y,A,alpha,a1,a2)

global_max_iter=100;   lasso_max_iter=100;
% Initialization and setting
M=size(y,1);   N=size(A,2);
X=zeros(N,1);  L=zeros(M,1);  Z=ones(M,1);
AT=A';  ATA =AT*A;
tau=eigs(ATA,1);  tauInv =1/tau;
beta=(20*M)/sum(sum(abs(y))); % betaInv=1./beta ;  % betaTauInv=betaInv*tauInv;
tolX=1e-6;  tolL=1e-6;  isConverged=0;  iter=0; 

while ~isConverged
    iter=iter+1;
    Xpre=X;    Lpre=L; 
   
    % Solve the low-rank regu;arized part
    tempM=(beta/(2+beta))*(A*X-y + (1/beta)*Z);
    [U,S,V]=svd(tempM,0);
    S=shrink(S,(alpha/(beta+2)));
    L=U*S*V';
   
   % Solve the sparse constraint part
   b=y-L+(1/beta).*Z;
   ATb=AT*b;
   X=FastLasso(ATA,ATb,X,tauInv,lasso_max_iter,a1,a2);

   % Update the Lagrange multiplier
   Z=Z+beta*(A*X-y-L);
   
   if (((norm(Xpre- X) < tolX*norm(Xpre)) && ...
       (norm(Lpre - L) < tolL*norm(Lpre))) || ...
       iter > global_max_iter) 
       isConverged = 1;
   end   

   iterationError=sum(sum(sqrt((A*X-y-L).^2)));
   fprintf('Error of iteration %d = %f \n', iter , iterationError);
end
end

