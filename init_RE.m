function [x, nu] = init_RE(dims)
rng(2011,'twister')
n=dims(1);
m=dims(2);
x=rand(n,1);
nu=rand(m,1);