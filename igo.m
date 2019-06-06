function y=igo(x,a,c)
[s1,s2]=size(x);
y=zeros(s1,s2);

if s1==1260
    nn=[42 30];
end

for i=1:s2
    t1=reshape(x(:,i),nn);
   fh=[-1 0 1;-2 0 2;-1 0 1]; fv=[-1 -2 -1;0 0 0;1 2 1]; % Sobel Filter can be change to others
   Gh=filter2(fh,t1); Gv=filter2(fv,t1);
   t2=tanh(a*(Gv./Gh-c));
    %t2=atan(a*(Gv./Gh-c));
   % t2= d*sigmf(Gv./Gh,[a c]);
   % t2 = softsign(Gv./Gh,[a c]);
    y(:,i)=reshape(t2,[s1,1]);
    y(find(isnan(y)==1))=0;
end
end