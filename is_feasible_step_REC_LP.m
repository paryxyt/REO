function [ flag ] = is_feasible_step_REC_LP(y,deltay,tbtl,varargin)
c=varargin{1};
K_params=varargin{2};
dims=varargin{3}; %[N; m; m1; n1; n2]

n1=dims(4)+dims(6); %variables in RER constraint
%n1=dims(3)+dims(4); %variables in RER constraint
n2=dims(5); %variables in orthant
m=dims(2); % equality constraints
m1=dims(3); % number of RER constraints

B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);
yp=y+tbtl*deltay;
Y1=yp(1:n1);
Y2=yp(n1+1:end);

if (min(B'*Y1>0) && min(C'*Y1>0)) && (max((B'*Y1).*log((B'*Y1)./(C'*Y1))+D'*Y1) <= 0) && (min(Y2) > 0)
    flag=1;
else
    flag=0;
end