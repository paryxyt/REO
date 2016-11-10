function [ flag ] = is_feasible_Phase1_REC_LP(y,varargin)
c=varargin{1};
K_params=varargin{2};
dims=varargin{3}; %[N; m; m1; n1; n2]

n1=dims(4)+dims(6); %variables in RER constraint
%n1=dims(3)+dims(4); %variables in RER constraint

n2=dims(5); %variables in orthant
m=dims(2); % equality constraints
m1=dims(3); % number of RER constraints

%the first n1-m1 components are the true variables, the last m1 are the slack variables 
B=K_params(1:n1-1,1:m1);
C=K_params(1:n1-1,m1+1:2*m1);
D=K_params(1:n1-1,2*m1+1:3*m1);
y_tr=y(1:n1-1);

% B=K_params(1:n1-m1,1:m1);
% C=K_params(1:n1-m1,m1+1:2*m1);
% D=K_params(1:n1-m1,2*m1+1:3*m1);
% y_tr=y(1:n1-m1);

z_tr=y(n1+1:end);

if (min(B'*(y_tr)>0) && min(C'*(y_tr)>0)) && (max((B'*y_tr).*log((B'*y_tr)./(C'*y_tr))+D'*y_tr) <= 0) && (min(z_tr)>0)
    flag=1;
else
    flag=0;
end

