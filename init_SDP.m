function [x, nu] = init_SDP(dims)
n=dims(1);
m=dims(2);
p=dims(3);

X=eye(p,p);
x=reshape(X,[p*p,1]);
nu=rand(m,1);