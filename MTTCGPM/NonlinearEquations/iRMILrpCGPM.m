% Matlab Model by Jianghua Yin (Jan.,2021, Nanning)
% Copyright (C) 2020 Jian Group
% All Rights Reserved
%
%% the derivative-free projection method (DFPM) for solving the 
%% unconstrained nonlinear monotone equations of the form
%   F(x)=0. 
%
function [Tcpu, NF, Itr, NormF] = iRMILrpCGPM(NO, method, model, x0,x0_old,para) 
 
format long

% start the clock
tic;

%% the number of itrations
% Itr=0;   

%% Initial information
[nprob,~] = init(NO);

%% the stopping criterion
epsilon = 1e-6;
epsilon1 = 1e-6;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;  
rho = para.rho;     % the compression ratio
fprintf('%s & %s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & alpha=%.4f & phi=%.4f  \n', ... 
    nprob,method,model,gamma,sigma,tau,rho);

 %% compute the search direction
Fk = feval(nprob,x0);   % evaluate the function value specified by nprob at x0在x0处计算nprob指定的函数值
NF = 1;
NormFxk = norm(Fk);
%x0_old = x0;
%x0_old1 = x0_old;
L1 = 0;

for k = 1:k_max
    
    if k==1 && NormFxk <= epsilon
        L1 = 1;
        NormF = NormFxk;    % the final norm of equations
        break;
    end
    
    % %% compute the inertial step %%
       if k==1
          eat=1;
       else   
          eat=1/k^2;
       end
     if x0 ~= x0_old        
     theta = min(0.01,eat/(norm(x0-x0_old)^2));
     else
      theta = 0.01;   
     end
    w0 = x0 + theta*(x0 - x0_old);
    Fwk = feval(nprob,w0);
    NF = NF + 1;
    NormFwk = norm(Fwk);   %求范数
    
    if NormFwk <= epsilon
        L1 = 1;
        NormF = NormFwk;   % the final norm of equations
        break;
    end
       
       %% compute the initial direction 计算初始方向%%
    if k==1
        dk = -Fwk;
    else
        
        % update the search direction
        switch method
             case 'iRMILrp-CGPM'
                yk=Fwk-Fwk_old;                          
                betak11=(Fwk'*yk)/(norm(dk)^2);                  
                betak2=(norm(Fwk)^2-(norm(Fwk))/(norm(Fwk))*(Fwk'*Fwk))/(norm(Fwk)^2+norm(Fwk)^2);
                betak1=max(0,betak11);                          
                if norm(yk)<= 0.1*norm(dk)
                     dk=-Fwk+betak1*dk ;
                else
                     dk=-1.5*Fwk+betak2*Fwk;
                end
                
        end
    end
  
    Normdk = norm(dk);
  
    if Normdk < epsilon1
        L1 = 1;
        NormF = NormFwk;
        break;
    end
  
    
    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk ≥ sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))*norm(dk)^2
    % model=3 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))/(1+norm(F(zk)))*norm(dk)^2
    % model=4 means -F(zk)'*dk ≥ sigma*tk*min(nu,norm(Fz_new,2))*norm(dk)^2
    if model==1
        t = gamma;
        z_new = w0 + t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        Normdk2 = norm(dk)^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0 + t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF + 1;
        end %%% End Armijo-type line search %%%
    elseif model==2
        t = 1;
        z_new = w0 + t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF+1;
        Normdk2 = norm(dk)^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*norm(Fz_new)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0 + t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF + 1;
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        z_new = w0 + t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF + 1;
        NormFzk = norm(Fz_new);
        Normdk2 = norm(dk)^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*NormFzk/(1+NormFzk)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0 + t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
            NormFzk = norm(Fz_new);
        end %%% End Armijo-type line search %%%
    elseif model==4   
        t = gamma;
        z_new = w0 + t*dk;
        Fz_new = feval(nprob,z_new);
        NF = NF + 1;
        NormFzk = norm(Fz_new);
        Normdk2 = norm(dk)^2;
        % check the Armijo-type line search condition
        while -Fz_new'*dk < sigma*t*max(0.01,min(0.8,NormFzk))*Normdk2/(1+Normdk2) && t>10^-6
            % the Armijo-type line search condition violated
            t = t*tau;
            z_new = w0 + t*dk;
            Fz_new = feval(nprob,z_new);
            NF = NF+1;
            NormFzk = norm(Fz_new);
        end %%% End Armijo-type line search %%%
    end 
    zk = z_new;
    Fzk = Fz_new;
    NormFzk = norm(Fzk);
    
    if NormFzk <= epsilon
        L1 = 1;
        NormF = NormFzk; % the final norm of equations
        break;
    end
    Fwz = Fzk'*(w0 - zk);
    xik = Fwz/NormFzk^2;
    % compute the next iteration 
    x1 = w0 - rho*xik*Fzk;
    
    Fk = feval(nprob,x1);
    NF = NF+1;
    NormFk = norm(Fk);
    
    if NormFk <= epsilon
        L1 = 1;
        NormF = NormFk;
        break;
    end

    % update the iteration
    dk_old = dk;
    Fwk_old = Fwk;
    x0_old = x0;
    x0 = x1;
    w0_old = w0;
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
