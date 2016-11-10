function [x, nu] = init_LP(dims)
n=dims(1);
m=dims(2);
x=rand(n,1);
nu=rand(m,1);