function [ deltaX deltaNu ] = infeasibleNewton_Phase1(Hf,A,b,gradf,X,nu,dims)
% Infeasible Newton System
n1=dims(1);
n2=dims(2);
m1=dims(3);
m2=dims(4);

M=[Hf A';
    A zeros(size(A,1),size(A,1))];
v=[gradf; A*X-b ];
%v=[gradf+A'*nu; A*X-b ];

delta=-pinv(M)*v;
deltaX=delta(1:n1+n2+m1);
w=delta(n1+n2+m1+1:end);

% Feasible Newton System
%v=[gradf; zeros(m2,1) ];
%delta=-M\v;
%deltaX=delta(1:n1+n2);

deltaNu=w-nu;
%deltaNu=delta(n2+1:end);

end

