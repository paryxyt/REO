function [F, gradF, H, fconstraint_orth] = barrier_LP(z,varargin)
c=varargin{1};
%% Define the barrier function

%%%%%%%%% Orthant part
fconstraint_orth=z;
F_orth=-sum(log(z)); %self-concordant barrier for z>=0
grad_orth=-(1./z);
H_orth=diag(1./(z.^2));

%%%%%%% Full barrier, gradient, Hessian
fconstraint=fconstraint_orth;
F = F_orth;
gradF =  grad_orth + c;

H =  H_orth;
