function [ X, nu, obj,exit_stat  ] = solver_primal( dims, A, b, c, Kparams, X_init, nu_init, barrier, is_feasible_step, varargin  )

if isempty(varargin) 
    phase1=0;
else
    if varargin{1}=='phase1'
        phase1=1;
        is_feasible_Phase1=varargin{2};
    else
        phase1=0;
    end
        
end
%%%%%%%% Phase 1 Needs to be properly handled %%%%%%%%%%%%%%%%%%

%% algorithm parameters
alpha=0.2;
%alpha=0.25;
beta=.5;
t=1; %initialization of t
mu=1.8;
eps=1e-6; %accuracy

newtonStep=@infeasibleNewton_Phase1; % specify the newton update formula here
residual=@PrimalOnlyNewtonResidual_Phase1; %specify how the residuals are computed
debug_history=0;

%input_data=@input_data4; % specify the input data here
%init=@init4; % specify the initialization of variables here
%barrier=@barrier_REC2; % specify the barrier function, its gradient and derivatives here
%is_feasible_step=@is_feasible_step_REC; %specify conditions for a step to be feasible (wrt domain)

%[dims,A,b,c]=input_data();
n=dims(1);
m=dims(2);
c0=c;


%% Initialization
X=X_init;
nu=nu_init;


X_hist=[];
fconstraint_hist=[];
s_hist=[];


%% Perform input checks


%% Outer iterations
r=ones(m+n,1);
outer_iter=0;
exit_stat=0;
while max(m/t,1/t) > eps
    outer_iter=outer_iter+1;
    fprintf('\n outer iteration %d \n',outer_iter);
    fprintf('-------------------- \n',outer_iter);
    
    %% Compute gradient and Hessian of t*objective + barrier function
    [F, gradf, Hf, fconstraint] = barrier(X,c,Kparams,dims);
     
       %% Newton Method
    iter_Newton=0;
    iter_NewtonMax=10000;
    
    while ((~isempty(A) & norm(A*X-b)>eps) | norm(r) > eps)
        iter_Newton = iter_Newton+1;
        if iter_Newton > iter_NewtonMax
            exit_stat=-1;
            X=[];
            obj=[];
            fprintf('\n Newton Method exceeded max iterations, problem seems infeasible. Exit status 1. \n');
            return;
        end
        
        [deltaX deltaNu] = newtonStep(Hf,A,b,gradf,X,nu,dims);
        

%%%%%%%%%%%%% Backtracking Line Search %%%%%%%%%%%
        tbtl=1;
        if is_feasible_step(X,deltaX,tbtl,c,Kparams,dims) %last three are optional args
            tbtl=1;
        else
            
            while ~(is_feasible_step(X,deltaX,tbtl,c,Kparams,dims))
                tbtl=beta*tbtl;
            end
        end
        if isnan(tbtl) | isinf(tbtl)
            fprintf('\n Warning: Backtracking Line Search halted. Moving to next Newton iteration \n')
            exit_stat=-1;
            tbtl=0;
            break;
        end
        iter=0;
        itermax_BTL=30;
        rold=residual(A,X,b,gradf,nu);
        X_p=X + tbtl*deltaX;
        nu_p=nu + tbtl*deltaNu;
        
        [~, gradf2, ~, ~] = barrier(X_p,c,Kparams,dims);
        r=residual(A,X_p,b,gradf2,nu_p);
        %r=residual(A,X_p,b,gradf+tbtl*Hf*deltaX,nu_p);
            
        while norm(r) > (1-alpha*tbtl)*norm(rold) 
            %pause
            iter=iter+1;
            tbtl=beta*tbtl;
            X_p = X + tbtl*deltaX;
            nu_p = nu + tbtl*deltaNu;
            if iter > itermax_BTL
                exit_stat=-1;
                fprintf('\n Backtracking line search failed. No Feasible Solution Found. Exit status 1. \n');
                X=[];
                obj=[];
                return
            end
            
             [~, gradf2, ~, ~] = barrier(X_p,c,Kparams,dims);
             r=residual(A,X_p,b,gradf2,nu_p);
        end
        %%%%%%%%%%%%%End Backtracking line search %%%%%%%%%%%%%%%%%%%%%%
        if rcond(Hf) < 1e-10 && phase1~=1
            break
        end

        X=X_p;
        nu=nu_p;
        [F, gradf, Hf, fconstraint] = barrier(X,c,Kparams,dims);
        r=residual(A,X,b,gradf,nu);
        
    fprintf('\n objective function value = %d',c0'*X);
    fprintf('\n equality constraint residual = %d \n',norm(A*X-b));
        
    
    % debug history ----------
    if debug_history==1
        X_hist =[X_hist X];
        fconstraint_hist=[fconstraint_hist fconstraint];
        s_hist=[s_hist X];
        if iter_Newton==3
            return
        end
    end
    if phase1==1 && norm(A*X-b) < eps && is_feasible_Phase1(X,c,Kparams,dims)
        obj=c0'*X;
        return
    end
    end % End of Newton step.
    t=mu*t;
    c=t*c;
    [F, gradf, Hf, fconstraint] = barrier(X,c,Kparams,dims);
    r=residual(A,X,b,gradf,nu);

    fprintf('\n Objective function value: %d \n',c0'*X)
    fprintf('\n fconstraint: %d \n',min(fconstraint))
    %pause
    obj=c0'*X;
end

end


