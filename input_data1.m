function [dims,A,b,c,Kparams]=input_data1()
% defines problem parameters
n=3;
m=1;
%A=randn(m,n);
A=[1 -1 1];
%b=rand(m,1);
c=rand(n,1);
b=zeros(m,1);
b=1;
dims=[n;m];
Kparams=[];