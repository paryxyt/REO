clear
clc
rng(2011,'twister')

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
 input_data=@input_data_RE; % specify the input data here
 init=@init_RE; % specify the initialization of variables here
 barrier=@barrier_REC2; % specify the barrier function, its gradient and derivatives here
 is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)

% Call the solver
[ y, obj, exit_status  ] = solver_primal( input_data, init, barrier, is_feasible_step  )