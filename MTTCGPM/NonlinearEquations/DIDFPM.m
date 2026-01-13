% Matlab Model by Jianghua Yin (Jan.,2021, Nanning)
% Copyright (C) 2020 Jian Group
% All Rights Reserved
%
%% the inertial-relaxed derivative-free projection method (IRDFPM) for solving the 
%% unconstrained nonlinear monotone equations of the form
%   F(x)=0. 
%
function [Tcpu,NF,Itr,NormF] =DIDFPM(NO,method,model,x0,para) 
 
format long

% start the clock
tic;

%% the number of itrations
% Itr=0;   

%% Initial information
[nprob,~]=init(NO);

%% the stopping criterion
epsilon=1e-6;
epsilon1=1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
beta = para.beta;% the coefficient of line search 
alpha1= para.alpha1;
thetak= para.thetak;
sigma = para.sigma;
tau = para.tau;         % the compression ratio
rho = para.rho;         % the relaxation factor 

fprintf('%s & %s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & alpha1=%.4f & rho=%.4f\n', ... 
    nprob,method,model,gamma,sigma,tau,alpha1,rho,thetak,beta);

%% compute the search direction
Fxk = feval(nprob,x0);   % evaluate the function value specified by nprob at x0
NF = 1;  
NormFxk = norm(Fxk);                   
x0_old = x0;
L1 = 0;
     
for k=1:k_max
    
    if k==1 && NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk; % the final norm of equations
        break; 
    end
    
    %% compute the inertial step %%
    u0 = x0 + beta*(x0 - x0_old);
    w0 = x0 + thetak*(x0 - x0_old);
    Fwk = feval(nprob,w0); 
    NF = NF+1;
    NormFwk = norm(Fwk);
    if NormFwk<=epsilon
        L1 = 1;
        NormF = NormFwk;   % the final norm of equations
        break; 
    end
    
    %% compute the initial direction %%
    if k==1
        dk = -Fwk;
    else
        % update the search direction
        switch method
            %case 'IRGPM'
                %dk = (-1)*Fwk; 
            case 'DIDFPM'
                sk = w0-w0_old;
                yk = Fwk-Fwk_old;
                gk = yk'*sk*yk+1e-3*(norm(sk))^2*sk; 
                theta_k = (sk'*sk)^2/(gk'*sk);
                dk = -theta_k*Fwk;
        end
    end
    Normdk = norm(dk);
    if Normdk<epsilon
        L1 = 1;
        NormF = NormFwk;
        break;
    end
    
    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk ¡Ý sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk ¡Ý sigma*tk*norm(F(zk))*norm(dk)^2
    % model=3 means -F(zk)'*dk ¡Ý sigma*tk*norm(F(zk))/(1+norm(F(zk)))*norm(dk)^2
    % model=4 means -F(zk)'*dk ¡Ý sigma*tk*min(nu,norm(Fz_new,2))*norm(dk)^2
    if model==1
        t = gamma;
        z_new = y0+t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        Normdk2 = Normdk^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = y0+t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
        end %%% End Armijo-type line search %%%
    elseif model==2
        t = gamma;
        z_new = w0+t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        Normdk2 = Normdk^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*norm(Fz_new)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0+t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        z_new = y0+t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        NormFzk = norm(Fz_new);
        Normdk2 = Normdk^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*NormFzk/(1+NormFzk)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = y0+t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
            NormFzk = norm(Fz_new);
        end %%% End Armijo-type line search %%%
    elseif model==4
        nu = 0.8;
        t = gamma;
        z_new = w0+t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        NormFzk = norm(Fz_new);
        Normdk2 = Normdk^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*min(nu,NormFzk)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0+t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
            NormFzk = norm(Fz_new);
        end %%% End Armijo-type line search %%%
    end 
    zk = z_new;
    Fzk = Fz_new;
    NormFzk = norm(Fzk);
    if NormFzk<=epsilon
        L1 = 1;
        NormF = NormFzk; % the final norm of equations
        break;
    end
    xik = Fzk'*(w0-zk)/NormFzk^2;
    % compute the next iteration 
    x1 = (1-alpha1)*u0 + alpha1*(w0-rho*xik*Fzk);
    Fxk = feval(nprob,x1);
    NF = NF+1;
    NormFxk = norm(Fxk);
    if NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk;
        break;
    end
    
    % update the iteration
    x0_old = x0;
    x0 = x1;
    w0_old = w0;
   Fwk_old = Fwk;
   % NormFwk_old =  NormFwk;
end
if L1==1
    Itr = k;
    Tcpu = toc;
else
    NF = NaN;
    Itr = NaN;
    Tcpu = NaN;
    NormF = NaN;
end
