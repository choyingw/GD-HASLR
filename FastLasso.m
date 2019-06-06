function [x]=FastLasso(AtA,Atb,x0,tauInv,lasso_max_iter,a1,a2)

tol_apg=1e-6;  isConverged=0;  iter=0;
x=x0;  t1=1; z=x ;

while ~isConverged

    iter=iter+1;
    xprev=x;
    
    temp1=z-tauInv*(AtA*z-Atb);
    beta=2./-pen_NIG(x,a1,a2);
    muTauInv=(1./beta).*tauInv;
    x=shrink(temp1,muTauInv);

    if (norm(xprev-x) < tol_apg*norm(xprev)) || iter>=lasso_max_iter
        isConverged=1;
    end
    
    t2=(1+sqrt(1+4*t1*t1))/2;
    z=x+((t1-1)/t2)*(x-xprev);
    t1=t2;
end

end

