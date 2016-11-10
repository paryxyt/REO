function [ flag ] = is_feasible_step_SDP(X,deltaX,tbtl,varargin)
p=round(sqrt(size(X,1)));
Z=reshape(X,[p,p]);
deltaZ=reshape(deltaX,[p,p]);
    if eig(Z+tbtl*deltaZ)>0
        flag=1;
    else
        flag=0;
    end