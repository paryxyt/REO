clear
clc

%% Linear Programming
 %input_data=@input_data_LP; % specify the input data here
 %init=@init_LP; % specify the initialization of variables here
 %barrier=@barrier_LP; % specify the barrier function, its gradient and derivatives here
 %is_feasible_step=@is_feasible_step_LP; %specify conditions for a step to be feasible (wrt domain)

%% SDP
%input_data=@input_data_SDP; % specify the input data here
%init=@init_SDP; % specify the initialization of variables here
%barrier=@barrier_SDP; % specify the barrier function, its gradient and derivatives here
%is_feasible_step=@is_feasible_step_SDP; %specify conditions for a step to be feasible (wrt domain)

%% Relative Entropy Optimization

 %input_data=@input_data_RE; % specify the input data here
 %init=@init_RE; % specify the initialization of variables here
 %barrier=@barrier_REC2; % specify the barrier function, its gradient and derivatives here
 %is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)

 %% Relative Entropy Optimization 
% Phase1
% Minimize:    sum(s)
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= s    (*)
%              A*y=b

% Phase2
% Minimize:    c'*y
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= 0    (*)
%              A*y=b

 %input_data=@input_data_RE_Phase1; % specify the input data here
 input_data=@input_data_RE; % specify the input data here
 [dims,A,b,c,Kparams]=input_data();
 n1=dims(1); % number of variables
 m=dims(2); % number of RE constraints
 m1=dims(3);  % number of equality constraints
 dims_phase1=[n1 m m1];
 A_phase1 = [A zeros(m,m1)];
 c_phase1=zeros(n1+m1,1);
 c_phase1(n1+1:n1+m1)=1; 
 
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
    s=2*(x1.*log(x1./x2)+x3);
    X_init=[X_init; s];
    nu_init=rand(m,1);
 
 Kparams_phase1=[Kparams; zeros(m1,3*m1)];
 Kparams_phase1(n1+1:n1+m1,2*m1+1:3*m1)=eye(m1,m1);
 %init=@init_RE_Phase1; % specify the initialization of variables here
 %[X_init, nu_init]=init(dims_phase1);
 barrier=@barrier_REC; % specify the barrier function, its gradient and derivatives here
 is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)
 is_feasible_Phase1=@is_feasible_Phase1_REC;

 
% Call the solver to solve Phase 1
[ y, nu, obj, exit_status  ] = solver_primal(dims_phase1,A_phase1,b,c_phase1,Kparams_phase1,X_init, nu_init, barrier, is_feasible_step,'phase1',is_feasible_Phase1)

% Call the solver to solve Phase 2
input_data=@input_data_RE; % specify the input data here
 [dims,A,b,c,Kparams]=input_data();
barrier=@barrier_REC; % specify the barrier function, its gradient and derivatives here
is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)
X_init=y(1:dims(1));
nu_init=nu;
clear y nu obj exit_status
[ y, nu, obj, exit_status  ] = solver_primal(dims,A,b,c,Kparams,X_init, nu_init, barrier, is_feasible_step)
%return
%% RER + LP
% Phase1
% Minimize:    sum(s)
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= s    (*)
%              z>=0   
%              A1*y + A2*z = b


% Phase2
% Minimize:    c'*y
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= 0    (*)
%              z>=0
%              A1*y + A2*z = b
%              A=[A1 A2]
clear
clc

input_data=@input_data_RE_LP; % specify the input data here
[dims,A,b,c,Kparams]=input_data();
 
 N=dims(1)+1+dims(2);
 %N=dims(1)+dims(4)+dims(2);
 n1=dims(1); % number of RE variables
 n2=dims(2); % LP variables
 m=dims(3); % number of equality constraints
 m1=dims(4);  % number of RE constraints
 dims_phase1=[N m m1 n1 n2 1];
 
 % initialization
    X_init=rand(n1,1);
    % here we assume that B, C are elementwise non-neg. If they are not, we
    % have to solve a linear program to find X_init that satisfied
    % B'*X_init > 0 and C'*X_init > 0
    z_init=rand(n2,1);
    B=Kparams(:,1:m1);
    C=Kparams(:,m1+1:2*m1);
    D=Kparams(:,2*m1+1:3*m1);
    x1=B'*X_init;
    x2=C'*X_init;
    x3=D'*X_init;
    
    s=max(-2*(x1.*log(x1./x2)+x3));
    %s=-2*(x1.*log(x1./x2)+x3);
    
    X_init=[X_init; s; z_init];
    nu_init=rand(m,1);

 A_phase1 = [A(:,1:n1) zeros(m,1) A(:,n1+1:n1+n2)];
 c_phase1=zeros(n1+n2+1,1);
 c_phase1(n1+1)=1; 
 Kparams_phase1=[Kparams; zeros(1,3*m1)];
 Kparams_phase1(n1+1,2*m1+1:3*m1)=-1;
 
 %A_phase1 = [A(:,1:n1) zeros(m,m1) A(:,n1+1:n1+n2)];
 %c_phase1=zeros(n1+n2+m1,1);
 %c_phase1(n1+1:n1+m1)=1; 
 %Kparams_phase1=[Kparams; zeros(m1,3*m1)];
 %Kparams_phase1(n1+1:n1+m1,2*m1+1:3*m1)=-eye(m1,m1);
  
 barrier=@barrier_REC_LP; % specify the barrier function, its gradient and derivatives here
 is_feasible_step=@is_feasible_step_REC_LP; %specify conditions for a step to be feasible (wrt domain)
 is_feasible_Phase1=@is_feasible_Phase1_REC_LP;

 %[F, gradF, H, fconstraint] = barrier_REC_LP(X_init,c_phase1,Kparams_phase1,dims_phase1)
 % Call the solver to solve Phase 1
[ y, nu, obj, exit_status  ] = solver_primal(dims_phase1,A_phase1,b,c_phase1,Kparams_phase1,X_init, nu_init, barrier, is_feasible_step,'phase1',is_feasible_Phase1)

% Call the solver to solve Phase 2
input_data=@input_data_RE_LP; % specify the input data here
 [dims,A,b,c,Kparams]=input_data();
barrier=@barrier_REC_LP; % specify the barrier function, its gradient and derivatives here
is_feasible_step=@is_feasible_step_REC_LP; %specify conditions for a step to be feasible (wrt domain)
X_init=[y(1:dims(1)); y(dims(1)+2:dims(1)+dims(2)+1)];
nu_init=nu;
clear y nu obj exit_status

N=dims(1)+dims(2);
 %N=dims(1)+dims(4)+dims(2);
 n1=dims(1); % number of RE variables
 n2=dims(2); % LP variables
 m=dims(3); % number of equality constraints
 m1=dims(4);  % number of RE constraints
 dims_phase2=[N m m1 n1 n2 0];
 

[ y, nu, obj, exit_status  ] = solver_primal(dims_phase2,A,b,c,Kparams,X_init, nu_init, barrier, is_feasible_step)

