clear
clc
rng(2007,'twister')

%Solver for REPs using Primal only method
% Minimize:    sum(s1)+sum(s2)
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= s1_i    (*)
%              A1*y + A2*z = b                         (**)
%              z + s2 >= 0
%              s1 >=0
%              s2 >=0
% Calls grad_f, PrimalOnlyNewtonResidual, infeasibleNewton



%% test example data
n1=1; % size of y
n2=1; % size of z
m1=1; % size of s1
m2=1; % num inequalities of form (**)
dims=[n1 n2 m1 m2];


% % there are m1 inequalities of the form (*)
% B=zeros(n1,m1);
% C=zeros(n1,m1);
% D=zeros(n1,m1);
% 
% B=[ 1 0 0 0 0 0; 
%     0 1 0 0 0 0;
%     0 0 1 0 0 0];
% B=B';
B=1;
% C=[1 .1 .1 0 0 0;
%     0 1  .2 0 0 0;
%     .1 .2 1 0 0 0 ];
% C=C';
C=1;
% D=[0 0 0   -.1    0   0;
%    0 0 0     0  -.1   0;
%    0 0 0     0    0 -.1];
% D=-D';
D=0;
% A1=rand(m2,n1); % there are m2 equalities of the form (**)
% A2=rand(m2,n2); % the dimension of z is n2
% b=zeros(m2,1);
% %b=rand(m2,1);
A1=1; A2=1; b=1;
A=[A1 A2 zeros(m2,m1) zeros(m2,n2)];
c=zeros(n1+n2+m1+n2,1);
c(n1+n2+1:n1+n2+m1+n2,1)=1;
% 
%% algorithm parameters
dims=[n1 n2 m1 m2];
%alpha=0.4;
alpha=0.25;
beta=.5;
t=1; %initialization of t
mu=5;
eps=1e-6; %accuracy

%% Initialization
y1=ones(n1,1);
z1=ones(n2,1);
s11=1*rand(m1,1);
s21=1*rand(n2,1);
X=[y1; z1; s11; s21];
nu=ones(m2,1);



%X=[ 5.2962    6.8909    7.6595 -101.8874 -437.1209 -491.6250  356.5116   43.9800    0.0000    0.0000    0.0000    0.0000 0 ]';
%% Perform input checks


%% Outer iterations
r=1;
outer_iter=0;
while m2/t > eps
    outer_iter=outer_iter+1;
    fprintf('\n outer iteration %d \n',outer_iter);
    fprintf('-------------------- \n',outer_iter);
    
    
    %% Compute gradient and Hessian of t*objective + barrier function
    [F, gradf, Hf, fconstraint] = inhomogeneous_barrier_Phase1(X,B,C,D,t,c,dims);
    
       %% Newton Method
    iter_Newton=0;
    iter_NewtonMax=99;
    while (norm(A*X-b)>eps | norm(r) > eps)
        iter_Newton = iter_Newton+1;
        if iter_Newton > iter_NewtonMax
            exit_stat=1;
            fprintf('\n Newton Method exceeded max iterations, problem seems infeasible. Exit status 1. \n');
            return;
        end
        [deltaX deltaNu] = infeasibleNewton_Phase1(Hf,A,b,gradf,X,nu,dims);
        
        %%%%%%%%%%%%% Backtracking Line Search
        % set tbtl
        y=X(1:n1);
        z=X(n1+1:n1+n2);
        s1=X(n1+n2+1:n1+n2+m1);
        s2=X(n1+n2+m1+1:n1+n2+m1+n2);
        deltay=deltaX(1:n1);
        deltaz=deltaX(n1+1:n1+n2);
        deltas1=deltaX(n1+n2+1:n1+n2+m1);
        deltas2=deltaX(n1+n2+m1+1:n1+n2+m1+n2);
        tbtl=1;
        %-By.*log(By./Cy)-Dy+s1
        
      %  if (min(B'*(y+deltay)>0) && min(C'*(y+deltay)>0) && min(s1+deltas1>0) && min(s2+deltas2>0) && min(z+s2+deltaz+deltas2>0))
      if (min(B'*(y+deltay)>0) && min(C'*(y+deltay)>0) && min(s1+deltas1>0) && min(s2+deltas2>0))
            tbtl=1;
        else
            %while (max(B'*(y+tbtl*deltay)<0) | max(C'*(y+tbtl*deltay)<0) | max(s1+tbtl*deltas1<0) | max(-B'*(y+tbtl*deltay)*(log((B'*(y+tbtl*deltay))/(C'*(y+tbtl*deltay)))) - D'*(y+tbtl*deltay)+s1+tbtl*deltas1 < 0) | max(s2+tbtl*deltas2<0) | max(z+s2+tbtl*deltaz+tbtl*deltas2<0))
            %while (max(B'*(y+tbtl*deltay)<0) | max(C'*(y+tbtl*deltay)<0) | max(s1+tbtl*deltas1<0)| max(s2+tbtl*deltas2<0) | max(z+s2+tbtl*deltaz+tbtl*deltas2<0))
            while (max(B'*(y+tbtl*deltay)<0) | max(C'*(y+tbtl*deltay)<0) | max(s1+tbtl*deltas1<0)| max(s2+tbtl*deltas2<0) ) 
            tbtl=beta*tbtl;
            end
        end
        if isnan(tbtl)|isinf(tbtl)
            fprintf('\n Warning: Backtracking Line Search halted. Moving to next Newton iteration \n')
            tbtl=0;
            break
        end
        iter=0;
        itermax_BTL=30;
        rold=PrimalOnlyNewtonResidual_Phase1(A,X,b,grad_f_Phase1(X,B,C,D,c,t,dims),nu);
        X_p=X + tbtl*deltaX;
        nu_p=nu + tbtl*deltaNu;
        r=PrimalOnlyNewtonResidual_Phase1(A,X_p,b,grad_f_Phase1(X_p,B,C,D,c,t,dims),nu_p);
        
        while norm(r) > (1-alpha*tbtl)*norm(rold);
            iter=iter+1;
            tbtl=beta*tbtl;
            X_p = X + tbtl*deltaX;
            nu_p = nu + tbtl*deltaNu;
            if iter > itermax_BTL
                exit_stat=1;
                fprintf('\n Backtracking line search failed. No Feasible Solution Found. Exit status 1. \n');
                return
            end
            r=PrimalOnlyNewtonResidual_Phase1(A,X_p,b,grad_f_Phase1(X_p,B,C,D,c,t,dims),nu_p);
        end
%%%%%%%%%%%%%End Backtracking line search %%%%%%%%%%%%%%%%%%%%%%

        X=X_p;
        z1=X;
        nu=nu_p;
        [F, gradf, Hf, fconstraint] = inhomogeneous_barrier_Phase1(X,B,C,D,t,c,dims);
        r=PrimalOnlyNewtonResidual_Phase1(A,X,b,gradf,nu);
        
%        fprintf('\t residual within Newton step: %d \n',norm(r))
    %X
    %deltas1
    %deltas2
    %pause
    end % End of Newton step.
    t=mu*t;
    r=PrimalOnlyNewtonResidual_Phase1(A,X,b,grad_f_Phase1(X,B,C,D,c,t,dims),nu);
    fprintf('\n Objective function value: %d \n',c'*X)
    
    %% If at any point feasible point produced, halt.
    
    %if min(norm(A*X-b)<eps) & min(fconstraint(1:m1)+X(n1+n2+1:n1+n2+m1)>0) & min(X(n1+1:n1+n2)>0)
    %    fprintf('\n strictly feasible solution found \n')
        %return
    %end

end


% X= [5.2962    6.8909    7.6595 -101.8874 -437.1209 -491.6250  356.5116
% 43.9800    0.0000    0.0000    0.0000    0.0000    0.0000] 
% is the true solution as determined using cvx_test_REP.m
