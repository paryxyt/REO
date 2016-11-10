function [dims,A,b,c,Kparams]=input_data2()
% defines problem parameters
p=3;
n=p*p;
m=1;
A=[];
for i=1:m
    A1=rand(p,p);
    A1=A1+A1';
    A=[A; transpose(vec(A1))];
end
C=rand(p,p)+10*eye(p,p);
C=C+C';
c=vec(C);
b=zeros(m,1);
b=1;
dims=[n;m;p];

Kparams=[];