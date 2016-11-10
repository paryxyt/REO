function [F, gradF, H, fconstraint] = barrier_REC(z,varargin)
% Inputs: c, K_params, dims are the objective function, cone parameters and dims
% respectively of the problem
c=varargin{1};
K_params=varargin{2};
dims=varargin{3};

n1=dims(1);
m=dims(2); %is zero
m1=dims(3);


B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);


%% Define the barrier function
y=z;
By=B'*y;
Cy=C'*y;
Dy=D'*y;

%%%%%%%%% EXPONENTIAL CONE PART %%%%%%%%%

fconstraint_exp=-By.*log(By./Cy)-Dy; %>=0
%Fcomp_exp =-log(By.*log(Cy./By)-Dy-s1)-log(By)-log(Cy); %old, mistake
Fcomp_exp =-log(-By.*log(By./Cy)-Dy)-log(By)-log(Cy);
F_exp=sum(Fcomp_exp); % the barrier function is sum of component barrier functions
gradFcomp=zeros(n1,m1);
Htcomp=zeros(n1,n1,m1);
for i=1:m1
    
    x2=By(i);
    x3=Cy(i);
    x1=Dy(i);
    L=[D(:,i) B(:,i) C(:,i)];
    
    
G=(x2*log(x2/x3)+x1);

gradFcomp(:,i) = L*((1/G)*[ -1; 
    -1-log(x2/x3);
    x2/x3]-[0; 1/x2; 1/x3]);

% Hessian
T1=-1/G;
T2=1+log(x2/x3);
T3=-(x2/x3);

Ht=zeros(3,3);
Ht(1,1)=1/(G)^2;
Ht(1,2)= 1/G^2*(log(x2/x3)+1);
Ht(1,3)=1/G^2*(-x2/x3);
Ht(2,1)=Ht(1,2);
Ht(2,2)= (1+log(x2/x3))/G^2 - (1/x2)*(1/G) + 1/x2^2;%T1/x2 + T2/G^2*(log(x2/x3)+1) + 1/x2^2;
Ht(2,3)=(1/x3)*(1/G)-(x2/x3)*(1+log(x2/x3))*(1/G^2);%T1*(-1/x3)+T2/G^2*(-x2/x3);
Ht(3,3)=(x2/x3)^2*(1/G)^2-x2/x3^2*(1/G)+1/x3^2;%T1*x2/x3^2+x2^2/x3^2*(1/G^2)+1/x3^2;
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