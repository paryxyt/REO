function [ r ] = PrimalOnlyNewtonResidual_Phase1( A,X,b,gradf,nu )
r_primal = A*X-b;
r_dual = gradf + A'*nu;
r=[r_primal; r_dual];
end

