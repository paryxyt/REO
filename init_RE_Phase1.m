function [x, nu] = init_RE_Phase1(dims)
n=dims(1);
m=dims(2);
x=rand(n,1);
x=x-10*[0 0 0 0 0 0 1 1 1]';
nu=rand(m,1);