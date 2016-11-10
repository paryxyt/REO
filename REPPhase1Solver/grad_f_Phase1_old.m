function [ gradf ] = grad_f_Phase1(X,B,C,D,c,t,dims)
n1=dims(1);
n2=dims(2);
m1=dims(3);
m2=dims(4);
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

end

grad_exp=[sum(gradFcomp,2);
    zeros(n2,1);
    gradFs1;
    zeros(m2,1)];


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

gradf = t*c + gradF; 
end

