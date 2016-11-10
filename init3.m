function [x, nu] = init3(dims)
n=dims(1);
m=dims(2);
x=ones(n,1);
nu=rand(m,1);