function y=igo3(x,a,c)
[s1,s2]=size(x);
y=zeros(s1,s2);
if s1==1260
    nn=[42 30];
end

for i=1:s2
    t1=reshape(x(:,i),nn);
    fh=[-1 0 1;-2 0 2;-1 0 1]; fv=[-1 -2 -1;0 0 0;1 2 1]; % Sobel Filter can be change to others
    Gh1=filter2(fh,t1); Gv1=filter2(fv,t1);
    Gh2=filter2(fh,Gh1); Gv2=filter2(fv,Gv1);
    Gh3=filter2(fh,Gh2); Gv3=filter2(fv,Gv2);
   % t2=atan(a*(Gv3./Gh3-c));
    t2=tanh(a*(Gv3./Gh3-c));
%    t2= d*sigmf(Gv3./Gh3,[a c]);
   % t2 = softsign(Gv3./Gh3,[a c]);
    y(:,i)=reshape(t2,[s1,1]);
    y(find(isnan(y)==1))=0;
end
end