function [ flag ] = is_feasible_step_REC(y,deltay,tbtl,varargin)
c=varargin{1};
K_params=varargin{2};
dims=varargin{3};

n1=dims(1);
m=dims(2); %is zero
m1=dims(3);


B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);
yp=y+tbtl*deltay;

if (min(B'*(y+tbtl*deltay)>0) && min(C'*(y+tbtl*deltay)>0)) && (max((B'*yp).*log((B'*yp)./(C'*yp))+D'*yp) <= 0)
    flag=1;
else
    flag=0;
end