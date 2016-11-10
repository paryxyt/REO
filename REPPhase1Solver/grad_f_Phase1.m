function [ gradF ] = grad_f_Phase1(X,B,C,D,c,t,dims)
n1 = dims(1); 
n2 = dims(2); 
m1 = dims(3); 
m2 = dims(4);
y=X(1:n1);
z=X(n1+1:n1+n2);
s1=X(n1+n2+1:n1+n2+m1);
%s2=X(n1+n2+m1+1:n1+n2+m1+n2);
By=B'*y;
Cy=C'*y;
Dy=D'*y;

%%%%%%%%% EXPONENTIAL CONE PART %%%%%%%%%

%fconstraint_exp=By.*log(By./Cy)-Dy-s1; % constraint function - is never used directly
fconstraint_exp=-By.*log(By./Cy)-Dy+s1; %>=0
%Fcomp_exp =-log(By.*log(Cy./By)-Dy-s1)-log(By)-log(Cy); %old, mistake
Fcomp_exp =-log(-By.*log(By./Cy)-Dy+s1)-log(By)-log(Cy);
%Fcomp_exp =-log(-By.*log(By./Cy)-Dy-s1)-log(By)-log(Cy);

F_exp=sum(Fcomp_exp); % the barrier function is sum of component barrier functions
gradFcomp=zeros(n1,m1);
for i=1:m1
    
    x2=By(i);
    x3=Cy(i);
    x1=Dy(i);
    L=[D(:,i) B(:,i) C(:,i)];
    s1i=s1(i);
    
G=(x2*log(x2/x3)+x1-s1i);

gradFcomp(:,i) = L*((1/G)*[ -1; 
    -1-log(x2/x3);
    x2/x3]-[0; 1/x2; 1/x3]);

gradFs1(i,1)=1/G;

% Hessian
end

grad_exp=[sum(gradFcomp,2);
    zeros(n2,1);
    gradFs1];



%%%%%%%%% Orthant part
fconstraint_orth=[z];
F_orth=-sum(log(z)); %self-concordant barrier for z>=0
%grad_orth=[zeros(n1,1); -(1./(z+s2));-(1./s1); -(1./(z+s2))-(1./s2)];
grad_orth=[zeros(n1,1); -(1./(z)); zeros(m1,1)];



%%%%%%% Full barrier, gradient, Hessian
fconstraint=[fconstraint_exp; fconstraint_orth];
F = F_exp + F_orth;
gradF = grad_exp + grad_orth + t*c;
%gradF = grad_exp + grad_orth;

end