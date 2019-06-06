
%% Main function of the method proposed in the  C. Y. Wu and J. J. Ding,¡§Occluded face recognition using low-rank regression with generalized gradient direction," under reviews of Pattern Recognition currently.
% For the data in AR database, pleas contact http://www2.ece.ohio-state.edu/~aleix/ARdatabase.html
% We split the whole database into neutral face image. scarf occlusion part, and sunglasses occlusion part from session 1 and session 2 respectively.
% After getting tha data, please resize the data image to 42x30 and construct the data as dictionaries to run the experiment in the paper

close all; clear all;
load('AR_1260_occu')
R=15; % How many people to extract as top residual 
r=1:300; %Specify the occlusion face number
accn=zeros(1,4); accnL1=zeros(1,4); accnL2=zeros(1,4); accnL3=zeros(1,4); %Size : n*i
u=7.3; v=0.51;% u: scale, v:shift
for n=1:1 % Controlling R index

trainDict= AR_1260_neutral; %trainDict= AR_1260_neutral12; 
train_faces=igo(trainDict,u,v);  train_facesL2=igo2(trainDict,u,v); train_facesL3=igo3(trainDict,u,v);
alpha = 100;  % balancing parameter 
people=100; % How many people
M=size(train_faces,1);
trainsize=size(train_faces,2); testsize=size(r,2);
factor_train=trainsize/people;  
count=0; countS1=0; countS2=0; countS3=0; 
a1=1;  a2=1e-6; % Parameters for NIG
    
for i=1:4 %Controlling the index of testing data 
    count=0; countS1=0; countS2=0; countS3=0;
    if i==1
        testDict=AR_1260_s1Scarf; 
    elseif i==2
        testDict=AR_1260_s1Sun;
    elseif i==3
        testDict=AR_1260_s2Scarf;
    else
        testDict=AR_1260_s2Sun;
    end
    
   factor_test=size(testDict,2)/people;
   
for l=1:testsize
  
 testImg=igo(testDict(:,r(l)),u,v);  testImgL2=igo2(testDict(:,r(l)),u,v); testImgL3=igo3(testDict(:,r(l)),u,v);

  disp('Solving Original');
  [X_recoveredL1,L_recoveredL1] = FastSolver(testImg,train_faces,alpha,a1,a2);
  disp('Solving Layer 2');
  [X_recoveredL2,L_recoveredL2] = FastSolver(testImgL2,train_facesL2,alpha,a1,a2);
  disp('Solving Layer 3');
  [X_recoveredL3,L_recoveredL3] = FastSolver(testImgL3,train_facesL3,alpha,a1,a2);
  
  residualL1=zeros(people,1); residualL2=zeros(people,1); residualL3=zeros(people,1);
  
  for j=1:people
       residualL1(j)=nuc_norm(testImg-train_faces(:,factor_train*(j-1)+1:factor_train*j)*X_recoveredL1(factor_train*(j-1)+1:factor_train*j)-L_recoveredL1);
       residualL2(j)=nuc_norm(testImgL2-train_facesL2(:,factor_train*(j-1)+1:factor_train*j)*X_recoveredL2(factor_train*(j-1)+1:factor_train*j)-L_recoveredL2);
       residualL3(j)=nuc_norm(testImgL3-train_facesL3(:,factor_train*(j-1)+1:factor_train*j)*X_recoveredL3(factor_train*(j-1)+1:factor_train*j)-L_recoveredL3);
  end
  [MinSL1,IdSL1]=min(residualL1); [MinSL2,IdSL2]=min(residualL2); [MinSL3,IdSL3]=min(residualL3);
  
  [MinL1,IdL1]=sort(residualL1,'ascend'); [MinL2,IdL2]=sort(residualL2,'ascend'); [MinL3,IdL3]=sort(residualL3,'ascend');
  I1=IdL1(1:R(n)).'; I2=IdL2(1:R(n)).'; I3=IdL3(1:R(n)).'; I=[I1 I2 I3];
  unqx=unique(I); his=histc(I,unqx); vcount=sort(his,'descend'); M=max(vcount);
  
  I_high=unqx(his==M);
  if length(I_high)>1 % Polling
      frequencyCount=zeros(length(I_high),1); b1=0; b2=0;
      for k=1:length(I_high)
        
          f1=find(IdL1==I_high(k)); f2=find(IdL2==I_high(k)); f3=find(IdL3==I_high(k));
          if isempty(f3)
              frequencyCount(k)= frequencyCount(k)+(R(n)+1);
          else
              frequencyCount(k)= frequencyCount(k)+f3; b1=1;
          end
          
          if isempty(f2)
              frequencyCount(k)= frequencyCount(k)+(R(n)+1);
          else
              frequencyCount(k)= frequencyCount(k)+f2; b2=1;
          end
          
          if b1&&b2
                  fr=0.3;
          else
                  fr=1;
          end
          
          if isempty(f1)
              frequencyCount(k)= frequencyCount(k)+fr*(R(n)+1);
          else
              frequencyCount(k)= frequencyCount(k)+fr*f1;
          end
          
      end
      [Mf,If]=min(frequencyCount);
      IdL=I_high(If);
  else
      IdL=I_high;
  end
  
  if(IdL==ceil(r(l)/factor_test))
        count=count+1; 
  end
      if(IdSL1==ceil(r(l)/factor_test))
          countS1=countS1+1; end
      if(IdSL2==ceil(r(l)/factor_test))
          countS2=countS2+1; end
      if(IdSL3==ceil(r(l)/factor_test))
          countS3=countS3+1; end
   
end
accn(n,i)=count/testsize; accnL1(n,i)=countS1/testsize;  accnL2(n,i)=countS2/testsize; accnL3(n,i)=countS3/testsize;
end
end
