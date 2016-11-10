function [ flag ] = is_feasible_step_LP(X,deltaX,tbtl,c,Kparams,dims)
    if X+tbtl*deltaX>0
        flag=1;
    else
        flag=0;
    end