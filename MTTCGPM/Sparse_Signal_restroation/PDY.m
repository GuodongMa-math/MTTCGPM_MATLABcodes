% Matlab Model by Jianghua Yin (July,2021, Nanning)
% Copyright (C) 2021 Jian Group
% All Rights Reserved
%
%% the inertial-relaxed derivative-free projection method (IRDFPM) for solving the 
%% unconstrained nonlinear monotone equations of the form
%   F(x)=0. 
%
function [x1,mses,Tcpu,NF,NormF] = PDY(A,b,lambda,true_x,method,model,para) 
 
format long

% if A is a matrix, we find out dimensions of y and x,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
  AT = @(x) (x'*A)'; %A'*x;
  A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% Precompute A'*b since it'll be used a lot
Atb = AT(b);    
x0 = Atb;
n = length(x0);

% initial point
xu0 =  x0.*(x0 >= 0);
xv0 = -x0.*(x0 <  0);
% from these two relations, x0 = u0-v0;

%% the stopping criterion
epsilon = 1e-6;
epsilon1 = 1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;         % the compression ratio
%alpha = para.alpha;     % the coefficient of the inertial step
rho = para.rho; % the relaxation factor 


fprintf('%s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & rho=%.4f\n', ... 
    method,model,gamma,sigma,tau,rho);

% start the clock
t0 = cputime;

% % define function handle
% Fu = @(x,xu) min(xu,AT(A(x))-Atb+lambda);
% Fv = @(x,xv) min(xv,-AT(A(x))+Atb+lambda);

%% compute the search direction
Ax0 = A(x0);
tempx = AT(Ax0)-Atb;
Fxu0 = min(xu0,tempx+lambda);
Fxv0 = min(xv0,-tempx+lambda);
NF = 1;  
NormFxk2 = Fxu0'*Fxu0+Fxv0'*Fxv0;    
NormFxk = sqrt(NormFxk2);
L1 = 0;
     
for k=1:k_max
    
    if k==1 && NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk; % the final norm of equations
        mses(k) = 1/n*norm(x0-true_x)^2;
        Tcpu(k) = cputime-t0;
        break; 
    end
    
    
    %% compute the initial direction %%
    if k==1
        dku = -Fxu0;
        dkv = -Fxv0;
    else
        % update the search direction
        switch method
            case 'DFSP1'
                dku = -Fxu0;
                dkv = -Fxv0; 
            case 'PDY'
             yku = Fxu0-Fxu0_old;
             ykv = Fxv0-Fxv0_old;
             Fxtdk = Fxu0'*dku+Fxv0'*dkv;
            dktyk = dku'*yku+dkv'*ykv;
            tku = 1+max(0,-dktyk/Normdk2);
            tkv = 1+max(0,-dktyk/Normdk2);
            wku = yku+tku*dku;
            wkv = ykv+tkv*dkv;
            dktwk = dku'*wku+dkv'*wkv;
%            betak = NormFxk2_old/dktwk;
           betak = NormFxk2/dktwk; 
            thetak =1- Fxtdk/dktwk;
            dku = (-1)*thetak*Fxu0+betak*dku; 
            dkv = (-1)*thetak*Fxv0+betak*dkv; 
 
        end
    end
    Normdk2 = dku'*dku+dkv'*dkv;
    Normdk = sqrt(Normdk2);
    if Normdk<epsilon1
        L1 = 1;
        NormF = NormFxk;
        mses(k) = 1/n*norm(yk-true_x)^2;
        Tcpu(k) = cputime-t0;
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
        zuk_new = xu0+t*dku;
        zvk_new = xv0+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = xuk+t*dku;
            zvk_new = xvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        end %%% End Armijo-type line search %%%
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
    elseif model==2
        t = gamma;
        zuk_new = xu0+t*dku;
        zvk_new = xv0+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*NormFzk_new*Normdk && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = xu0+t*dku;
            zvk_new = xv0+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        zuk_new = xuk+t*dku;
        zvk_new = xvk+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*NormFzk_new/(1+NormFzk_new)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = xuk+t*dku;
            zvk_new = xvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
        end %%% End Armijo-type line search %%%
    else
        nu = 0.8;
        t = gamma;
        zuk_new = xuk+t*dku;
        zvk_new = xvk+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*min(nu,NormFzk_new)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = xuk+t*dku;
            zvk_new = xvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
        end %%% End Armijo-type line search %%%
    end
    xu0_old = xu0;
    xv0_old = xv0;
    zuk = zuk_new;
    zvk = zvk_new;
    zk = zk_new;
    Fzuk = Fzuk_new;
    Fzvk = Fzvk_new;
    Fxu0_old = Fxu0;
    Fxv0_old = Fxv0;
    NormFxk_old = NormFxk;       
    NormFzk2 = NormFzk_new2;
    NormFzk = NormFzk_new;
    if NormFzk<=epsilon
        L1 = 1;
        NormF = NormFzk; % the final norm of equations
        mses(k) = 1/n*norm(zk-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    Fzktykzk = Fzuk'*(xu0-zuk)+Fzvk'*(xv0-zvk);
    xik = Fzktykzk/NormFzk2;
    % compute the next iteration 
    xu1 = xu0-rho*xik*Fzuk;
    xv1 = xv0-rho*xik*Fzvk;
    x1 = xu1-xv1;
    mses(k) = 1/n*norm(x1-true_x)^2;
    Ax1 = A(x1);
    tempx = AT(Ax1)-Atb;
    Fxu0 = min(xu1,tempx+lambda);
    Fxv0 = min(xv1,-tempx+lambda);
    NF = NF+1;
    NormFx2 = Fxu0'*Fxu0+Fxv0'*Fxv0;
    NormFxk = sqrt(NormFx2);
    if NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk;
        mses(k) = 1/n*norm(x1-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    
    % update the iteration
    xu0 = xu1;
    xv0 = xv1;
    Tcpu(k) = cputime-t0;
    disp(NF);
end
if L1~=1
    NF = NaN;
    NormF = NaN;
end
