function [y, nu, obj] = call_solver_REO(dims,A,b,c,Kparams)
% Relative Entropy Optimization 
% Phase1
% Minimize:    s
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= s    (*)
%              A*y=b
%
% Phase2
% Minimize:    c'*y
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= 0    (*)
%              A*y=b

 n1=dims(1); % number of variables
 m=dims(2); % number of equality constraints
 m1=dims(3);  % number of RE constraints
 ns=1;
 N=n1+1;
 dims_phase1=[N m m1 n1];
 A_phase1 = [A zeros(m,ns)];
 c_phase1=zeros(n1+ns,1);
 c_phase1(n1+1:n1+ns)=1; 
 
 % initialization
    X_init=rand(n1,1);
    % here we assume that B, C are elementwise non-neg. If they are not, we
    % have to solve a linear program to find X_init that satisfied
    % B'*X_init > 0 and C'*X_init > 0
    B=Kparams(:,1:m1);
    C=Kparams(:,m1+1:2*m1);
    D=Kparams(:,2*m1+1:3*m1);
    x1=B'*X_init;
    x2=C'*X_init;
    x3=D'*X_init;
    s=max(-2*(x1.*log(x1./x2)+x3));
    X_init=[X_init; s];
    nu_init=rand(m,1);
 
 Kparams_phase1=[Kparams; zeros(ns,3*m1)];
 Kparams_phase1(n1+1:n1+ns,2*m1+1:3*m1)=-1;
 %init=@init_RE_Phase1; % specify the initialization of variables here
 %[X_init, nu_init]=init(dims_phase1);
 barrier=@barrier_REC; % specify the barrier function, its gradient and derivatives here
 is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)
 is_feasible_Phase1=@is_feasible_Phase1_REC;

 
% Call the solver to solve Phase 1
[ y, nu, obj, exit_status  ] = solver_primal(dims_phase1,A_phase1,b,c_phase1,Kparams_phase1,X_init, nu_init, barrier, is_feasible_step,'phase1',is_feasible_Phase1)
y0=y;
%return

% Call the solver to solve Phase 2
%input_data=@input_data_RE; % specify the input data here
% [dims,A,b,c,Kparams]=input_data();
barrier=@barrier_REC; % specify the barrier function, its gradient and derivatives here
is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)
X_init=y(1:dims(1));
nu_init=nu;
clear y nu obj exit_status
[ y, nu, obj, exit_status  ] = solver_primal(dims,A,b,c,Kparams,X_init, nu_init, barrier, is_feasible_step)

