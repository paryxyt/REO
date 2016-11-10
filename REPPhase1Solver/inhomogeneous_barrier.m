function [F, gradF, H, fconstraint_exp, fconstraint_orth] = inhomogeneous_barrier(X,B,C,D,dims)
%% Define the barrier function
% There are TWO parts to the barrier function
% rel ent part and orthant part
% test points
%B=[1 1.1; 0 0.1 ; 0 0]; C=[0 0; 1 .1; 0 .05]; D=[0 0; 0 0; 1 1]; 
%y=[1.2; 1.7; -.3]; z=[1; 1]; tau=1;
[n1 n2 m1 m2] = dims;
y=X(1:n1);
z=X(n1+1:n1+n2);
s1=X(n1+n2+1:n1+n2+m1);
s2=X(n1+n2+m1+1:n1+n2+m1+n2);
By=B'*y;
Cy=C'*y;
Dy=D'*y;

%%%%%%%%% EXPONENTIAL CONE PART %%%%%%%%%

fconstraint_exp=By.*log(By./Cy)-Dy-s1; % constraint function - is never used directly
Fcomp_exp =-log(By.*log(Cy./By)-Dy-s1)-log(By)-log(Cy);
F_exp=sum(Fcomp_exp); % the barrier function is sum of component barrier functions
gradFcomp=zeros(n1,m1);
Htcomp=zeros(n1,n1,m1);
for i=1:m1
    
    Byi=By(i);
    Cyi=Cy(i);
    Dyi=Dy(i);
    L=[D(:,i) C(:,i) B(:,i)];
    s1i=s1(i);
    
gradFcomp(:,i) = L*[ 1/(Byi*log(Cyi/Byi)-Dyi-s1i);
    -1/(Byi*log(Cyi/Byi)-Dyi-s1i)*(1/Cyi)-(1/Cyi);
    (1-log(Cyi/Byi))/(Byi*log(Cyi/Byi)-Dyi-s1i)-(1/Byi)
];

gradFs1(i,1)=1/(Byi*log(Cyi/Byi)-Dyi-s1i);


% Hessian
x1=Dyi+s1i;
x2=Cyi;
x3=Byi;
T1=x3*log(x2/x3)-x1;
T2=1-log(x2/x3);

Ht=zeros(3,3);
Ht(1,1)=1/(T1)^2;
Ht(1,2)= - (x3/x2)*(1/T1^2);
Ht(1,3)=T2/(T1^2);
Ht(2,1)=Ht(1,2);
Ht(2,2)=1/(x2)^2*(1+(1/T1)) + (1/x2)*(1/T1^2)*(x3/x2);
Ht(2,3)=-(1/x2)*Ht(1,3);
Ht(3,1)=Ht(1,3);
Ht(3,2)=Ht(2,3);
Ht(3,3)=1/x3^2 + (1/T1^2)*((1/x3)*T1 + T2^2);
Htcomp(:,:,i)=L*Ht*L';

Hzs1(:,:,i)=L*[ 1/(Byi*log(Cyi/Byi)-Dyi-s1i).^2;
    (-1/(Byi*log(Cyi/Byi)-Dyi-s1i)^2)*(1/Cyi);
    (1-log(Cyi/Byi))/(Byi*log(Cyi/Byi)-Dyi-s1i)^2
];


Hs1s1(:,:,i)=1/(Byi*log(Cyi/Byi)-Dyi-s1i)^2;
end

grad_exp=[sum(gradFcomp,2);
    zeros(n2,1);
    gradFs1;
    zeros(m2,1)];


H_e1=sum(Htcomp,3);
Hzs1_full=sum(Hzs1,3);
Hs1z_full=transpose(sum(Hzs1,3));
Hs1s1_full=sum(Hs1s1,3);

H_exp=zeros(n1+n2+m1+n2,n1+n2+m1+n2); 

% Derivatives order in H
%[yy  yz  ys1  ys2
% zy  zz  zs1  zs2
% s1y s1z s1s1 s1s2
% s2y s2z s2s1 s2s2];

H_exp(1:n1,1:n1)=H_e1;
H_exp(n1+n2+1:n1+n2+m1,1:n1)=Hzs1_full;
H_exp(n1+n2+1:n1+n2+m1,n1+n2+1:n1+n2+m1)=Hs1z_full;
H_exp(1:n1,n1+n2+1:n1+n2+m1)=Hs1s1_full;
%%%%%%%%% Orthant part
fconstraint_orth=[z+s2; s1; s2];
F_orth=-sum(log(z+s2))-sum(log(s1))-sum(log(s2)); %self-concordant barrier for z>=0
grad_orth=[zeros(n1,1); -(1./(z+s2));-(1./s1); -(1./(z+s2))-(1./s2)];

H_orth=zeros(n1+n2+m1+n2,n1+n2+m1+n2); 

H_orth(n1+1:n1+n2,n1+1:n1+n2)= diag(1./(z+s2).^2);%zz

H_orth(n1+1:n1+n2,n1+n2+m1+1:n1+n2+m1+n2)= diag(1./(z+s2).^2); %zs2
H_orth(n1+n2+m1+1:n1+n2+m1+n2,n1+1:n1+n2)= diag(1./(z+s2).^2);%s2z

H_orth(n1+n2+1:n1+n2+m1,n1+n2+1:n1+n2+m1)= diag(1./s1.^2);%s1s1

H_orth(n1+n2+m1+1:n1+n2+m1+m2,n1+n2+m1+1:n1+n2+m1+m2)= diag(1./(z+s2).^2)+diag(1./s2.^2);%s2s2


%%%%%%% Full barrier, gradient, Hessian
fconstraint=[fconstraint_exp; fconstraint_orth];
F = F_exp + F_orth;
gradF = grad_exp + grad_orth;
H = H_exp + H_orth;  
