function [ deltaX deltaNu ] = infeasibleNewton_Phase1(Hf,A,b,gradf,X,nu,dims)
% Infeasible Newton System
n=dims(1);
m=dims(2);

M=[Hf A';
    A zeros(size(A,1),size(A,1))];
if isempty(A)
    q1=zeros(size(gradf,1),1);
    q2=[];
else
    q1=A'*nu;
    q2=A*X-b;
end

v=[gradf + q1; q2 ];
%v=[gradf+A'*nu; A*X-b ];

delta=-pinv(M)*v;
deltaX=delta(1:n);
deltaNu=delta(n+1:end);
%w=delta(n+1:end);

%deltaNu=w-nu;
%deltaNu=delta(n2+1:end);

end

