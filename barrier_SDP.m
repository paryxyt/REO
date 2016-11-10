function [F, gradF, H, fconstraint] = barrier_SDP(z,varargin)
c=varargin{1};

%% Define the barrier function
p=round(sqrt(size(z,1)));
Z=reshape(z,[p,p]);
%%%%%%%%% SDP part
fconstraint_sdp=svd(Z);
F_sdp=-log(det(Z)); %self-concordant barrier for z>=0
grad_sdp=-vec(inv(Z));
H_sdp=kron(inv(Z),inv(Z));

%%%%%%% Full barrier, gradient, Hessian
fconstraint=fconstraint_sdp;
F = F_sdp;
gradF =  grad_sdp + c;
H =  H_sdp;
