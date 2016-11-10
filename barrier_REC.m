function [F, gradF, H, fconstraint] = barrier_REC(Y,varargin)
% Inputs: c, K_params, dims are the objective function, cone parameters and dims
% respectively of the problem
c=varargin{1};
K_params=varargin{2};
dims=varargin{3};

n1=dims(1); %number of variables
m=dims(2); % number of equality constraints
m1=dims(3); % number of RE constraints


B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);


%% Define the barrier function

By=B'*Y;
Cy=C'*Y;
Dy=D'*Y;

%%%%%%%%% EXPONENTIAL CONE PART %%%%%%%%%

fconstraint_exp=-By.*log(By./Cy)-Dy; %>=0
Fcomp_exp =-log(-By.*log(By./Cy)-Dy)-log(By)-log(Cy);
F_exp=sum(Fcomp_exp); % the barrier function is sum of component barrier functions
gradFcomp=zeros(n1,m1);
Htcomp=zeros(n1,n1,m1);
for i=1:m1
    
    x=By(i);
    y=Cy(i);
    z=Dy(i);
    L=[B(:,i) C(:,i) D(:,i)];
    
   

gradFcomp(:,i) = L*[(-2*x*log(x/y)-x-z)/(x^2*log(x/y)+x*z);
    (1/y)*(x/(x*log(x/y)+z)-1);
    -1/(x*log(x/y)+z)];

% Hessian

Ht=zeros(3,3);
Ht(1,1)=(2*x^2*(log(x/y))^2 + x^2 + x*(x+2*z)*log(x/y) -x*z + z^2 )/(x*(x*log(x/y)+z))^2;
Ht(1,2)= (z-x)/(y*(x*log(x/y)+z)^2);
Ht(1,3)=(x-z)/(x*(x*log(x/y)+z)^2) + 1/(x*(x*log(x/y)+z));
Ht(2,1)=Ht(1,2);
Ht(2,2)= x^2/(y*(x*log(x/y)+z))^2 - x/(y^2*(x*log(x/y)+z)) + 1/y^2;
Ht(2,3)=-x/(y*(x*log(x/y)+z)^2);
Ht(3,3)=1/(x*log(x/y)+z)^2;
Ht(3,1)=Ht(1,3);
Ht(3,2)=Ht(2,3);
Htcomp(:,:,i)=L*Ht*L';
end

grad_exp=sum(gradFcomp,2);
H_e1=sum(Htcomp,3);
H_exp(1:n1,1:n1)=H_e1;


%%%%%%% Full barrier, gradient, Hessian
fconstraint=[fconstraint_exp];
F = F_exp;
gradF = grad_exp + c;
H = H_exp;  


end