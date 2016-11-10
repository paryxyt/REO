function [ r ] = PrimalOnlyNewtonResidual_Phase1( A,X,b,gradf,nu )
if isempty(A)
    r_primal=[];
    q2=zeros(size(gradf,1),1);
else
    r_primal = A*X-b;
    q2=A'*nu;
end
r_dual = gradf + q2;
r=[r_dual; r_primal];
end

