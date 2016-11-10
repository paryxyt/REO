function [ flag ] = is_feasible_Phase1_REC(y,varargin)
c=varargin{1};
K_params=varargin{2};
dims=varargin{3};

n1=dims(1); %number of variables
m=dims(2); %number of equality constraints
m1=dims(3); %number of RE constraints

%the first n1-m1 components are the true variables, the last m1 are the slack variables 
B=K_params(1:n1-1,1:m1);
C=K_params(1:n1-1,m1+1:2*m1);
D=K_params(1:n1-1,2*m1+1:3*m1);
y_tr=y(1:n1-1);

max((B'*y_tr).*log((B'*y_tr)./(C'*y_tr))+D'*y_tr)

if (min(B'*(y_tr)>0) && min(C'*(y_tr)>0)) && (max((B'*y_tr).*log((B'*y_tr)./(C'*y_tr))+D'*y_tr) <= 0)
    flag=1;
else
    flag=0;
end