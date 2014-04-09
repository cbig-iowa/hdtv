function mask = randSampling(n1,n2,n3,tim)
%tim--undersampling factor

len=n1*n2*n3;
pts=round(tim*len);

Stemp=zeros(len,1);
R=randperm(len);
U=ones(len,1);
Stemp(R(1:pts))=U(R(1:pts));
mask=reshape(Stemp,n1,n2,n3);
