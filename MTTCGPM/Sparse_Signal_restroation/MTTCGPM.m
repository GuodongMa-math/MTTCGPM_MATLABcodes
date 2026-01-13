% Matlab Model by Jianghua Yin (July,2021, Nanning)
% Copyright (C) 2021 Jian Group
% All Rights Reserved
%
%% the inertial-relaxed derivative-free projection method (IRDFPM) for solving the 
%% unconstrained nonlinear monotone equations of the form
%   F(x)=0. 
%
function [x1,mses,Tcpu,NF,NormF] = MTTCGPM(A,b,lambda,true_x,method,model,para) 
 
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

xu0 =  x0.*(x0 >= 0);
xv0 = -x0.*(x0 <  0);

xu0_old =  x0.*(x0 >= 0);
xv0_old = -x0.*(x0 <  0);

xu0_old1 =  x0.*(x0 >= 0);
xv0_old1 = -x0.*(x0 <  0);

%% the stopping criterion
epsilon = 1e-6;
epsilon1 = 1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;         % the compression ratio
rho = para.rho;         % the relaxation factor ËÉ³ÚÒò×Ó
mu = para.mu; 
xi = para.xi;  
eat = para.eat; 

fprintf('%s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & rho=%.4f  & mu=%.4f & xi=%.4f & eat=%.4f\n', ... 
    method,model,gamma,sigma,tau,rho,mu,xi,eat);

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
    
    %% compute the inertial step %%
    Normxkxk12 =(xu0-xu0_old)'*(xu0-xu0_old)+(xv0-xv0_old)'*(xv0-xv0_old);
    Normxkxk1 = sqrt(Normxkxk12 );
    ak=min(0.1,1/(k^2*Normxkxk1));
    Normxk1xk22 =(xu0_old-xu0_old1)'*(xu0_old-xu0_old1)+(xv0_old-xv0_old1)'*(xv0_old-xv0_old1);
    Normxkx1k2= sqrt(Normxk1xk22 );
    bk=min(0.01,1/(k^4*Normxkx1k2));
    yuk = xu0+ak*(xu0-xu0_old)+bk*(xu0_old-xu0_old1);
    yvk = xv0+ak*(xv0-xv0_old)+bk*(xv0_old-xv0_old1);
    yk = yuk-yvk;
    Ayk = A(yk);
    tempy = AT(Ayk)-Atb;
    Fyuk = min(yuk,tempy+lambda);
    Fyvk = min(yvk,-tempy+lambda);
    NF = NF+1;
    NormFyk2 = Fyuk'*Fyuk+Fyvk'*Fyvk;
    NormFyk = sqrt(NormFyk2);
    if NormFyk<=epsilon
        L1 = 1;
        NormF = NormFyk;   % the final norm of equations
        mses(k) = 1/n*norm(yk-true_x)^2;
        Tcpu(k) = cputime-t0;
        break; 
    end
    
    %% compute the initial direction %%
    if k==1
        dku = -Fyuk;
        dkv = -Fyvk;
    else
        % update the search direction
        switch method
            case 'ALgorithm2'
                dku = -Fyuk;
                dkv = -Fyvk; 
            case 'MTT-CGPM'
             uuk = Fyuk - Fyuk_old ;
             uvk = Fyvk - Fyvk_old ;
             puk= Fyuk_old;
             pvk= Fyvk_old;
             Normuk2 = uuk'*uuk+uvk'*uvk;
             Normuk = sqrt(Normuk2);
              Normpk2 = puk'*puk+pvk'*pvk;
             Normpk = sqrt(Normpk2);
             tuk=Normuk/Normdk-min(0,dku'*Fyuk/Normdk2);
             tvk=Normuk/Normdk-min(0,dkv'*Fyvk/Normdk2);
             t= tuk-tvk;
             wku1=uuk+tuk*dku;
             wkv1=uvk+tvk*dkv; 
             Normwk12 = wku1'*wku1+wkv1'*wkv1;
             Normwk1=sqrt(Normwk12);
             suk=tuk'*uuk;
             svk=tvk'*uvk;
             Normsk12 = suk'*suk+svk'*svk;
             Normsk1=sqrt(Normsk12);
             deteku=Fyuk'*wku1/Normwk12;
             detekv=Fyvk'*wkv1/Normwk12;
             betak11u=(NormFyk2-(NormFyk/Normpk)*Fyuk'*puk)/max(mu*Normdk*Normwk1+NormFyk2_old,dku'*wku1);
             betak11v=(NormFyk2-(NormFyk/Normpk)*Fyvk'*uvk)/max(mu*Normdk*Normwk1+NormFyk2_old,dkv'*wkv1);             
             betak2u=Fyuk'*puk/(NormFyk2+Normpk2)+(Normwk1/Normsk1)*suk'*Fyuk/(NormFyk*Normpk);
             betak2v=Fyvk'*pvk/(NormFyk2+Normpk2)+(Normwk1/Normsk1)*svk'*Fyvk/(NormFyk*Normpk);
             betak1u=min(0,betak11u);
             betak1v=min(0,betak11v);            
             theatu1=max(10,min(10^30,1+betak2u*Fyuk'*puk/NormFyk2));
             theatv1=max(10,min(10^30,1+betak2v*Fyvk'*pvk/NormFyk2));           
                if k>=2 && Normwk1>= xi*NormFyk 
                     dku=-Fyuk+betak1u*dku+(2-eat)*deteku*wku1 ;
                     dkv=-Fyvk+betak1v*dkv+(2-eat)*detekv*wkv1 ;
                else
                     dku=-theatu1*Fyuk+betak2u*puk;
                     dkv=-theatv1*Fyvk+betak2v*pvk;
                end
        end
    end
    Normdk2 = dku'*dku+dkv'*dkv;
    Normdk = sqrt(Normdk2);
    if Normdk<epsilon1
        L1 = 1;
        NormF = NormFyk;
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
            zuk_new = xu0+t*dku;
            zvk_new = xv0+t*dkv;
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
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*NormFzk_new*Normdk2 && t>10^-6  
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
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        zuk_new = xu0+t*dku;
        zvk_new = xv0+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*NormFzk_new/(1+NormFzk_new)*Normdk2 && t>10^-6  
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
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        end %%% End Armijo-type line search %%%
    else
         t = gamma;
        zuk_new = yuk+t*dku;
        zvk_new = yvk+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        eta_k=max(0.0001,min(0.8, NormFzk_new));
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*eta_k*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = yuk+t*dku;
            zvk_new = yvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
             eta_k=max(0.0001,min(0.8, NormFzk_new));
        end %%% End Armijo-type line search %%%     
    end
    zuk = zuk_new;
    zvk = zvk_new;
    zk = zk_new;
    Fzuk = Fzuk_new;
    Fzvk = Fzvk_new;
    NormFzk2 = NormFzk_new2;
    NormFzk = NormFzk_new;
    if NormFzk<=epsilon
        L1 = 1;
        NormF = NormFzk; % the final norm of equations
        mses(k) = 1/n*norm(zk-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    Fzktykzk = Fzuk'*(yuk-zuk)+Fzvk'*(yvk-zvk);
    xik = Fzktykzk/NormFzk2;
    % compute the next iteration 
    xu1 = yuk-rho*xik*Fzuk;
    xv1 = yvk-rho*xik*Fzvk;
    x1 = xu1-xv1;
    mses(k) = 1/n*norm(x1-true_x)^2;
    Ax1 = A(x1);
    tempx = AT(Ax1)-Atb;
    Fxu = min(xu1,tempx+lambda);
    Fxv = min(xv1,-tempx+lambda);
    NF = NF+1;
    NormFx2 = Fxu'*Fxu+Fxv'*Fxv;
    NormFxk = sqrt(NormFx2);
    if NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk;
        mses(k) = 1/n*norm(x1-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    
    % update the iteration
    xu0_old1 = xu0_old;
    xu0_old = xu0;
    xu0 = xu1;
    xv0_old1 = xv0_old;
    xv0_old = xv0;
    xv0 = xv1;
    Fyuk_old = Fyuk;
    Fyvk_old = Fyvk;
    NormFyk2_old = NormFyk2;
    Tcpu(k) = cputime-t0;
end
if L1~=1
    NF = NaN;
    NormF = NaN;
end
