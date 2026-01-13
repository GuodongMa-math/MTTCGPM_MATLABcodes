
clear all
close all
rng(2016);
ITR_max = 2000;

% parameters for TSIDFPM
para1.Itr_max = ITR_max;
para1.gamma =1;         % the initial guess
para1.sigma = 0.001;      % the coefficient of line search 
para1.tau =0.3;         % the compression ratio
para1.rho = 1;           % the relaxation factor 
     
% parameters for DIDFPM
para2.Itr_max = ITR_max;
para2.gamma = 0.3;
para2.sigma = 0.001;      % the coefficient of line search 
para2.tau = 0.3;         % the compression ratio
para2.alpha = 0.28;       % the coefficient of inertial step
para2.rho = 1.8;      % the relaxation factor 

% parameters for MTT-CGPM
para3.Itr_max = ITR_max;
para3.gamma = 0.45;         % the initial guess
para3.sigma = 0.0001;    % the c oefficient of line search
para3.tau = 0.4;         % the compression ratio
para3.rho = 1.9;
para3.mu = 3;
para3.xi= 1;
para3.eat= 2.07;

% parameters for iRMILrp-CGPM
para4.Itr_max = ITR_max;
para4.gamma = 0.6;         % the initial guess
para4.sigma = 0.01;      % the coefficient of line search
para4.tau = 0.6;         % the compression ratio
para4.rho = 1.1;    

% parameters for PDY
para5.Itr_max = ITR_max;
para5.gamma = 1;         % the initial guess
para5.sigma = 0.01;      % the coefficient of line search
para5.tau = 0.5;         % the compression ratio
para5.rho = 1;

% n is the original signal length
n = 2^12;

% k is number of observations to make
m = 2^10;

% number of spikes to put down
n_spikes = 64;

% random +/- 1 signal
f = zeros(n,1);
q = randperm(n);
% f(q(1:n_spikes)) = sign(randn(n_spikes,1));
f(q(1:n_spikes)) = randn(n_spikes,1);

% measurement matrix
disp('Creating measurement matrix...');
R = randn(m,n);

% orthonormalize rows
R = orth(R')';

if n == 8192  
   % in this case, we load a precomputed
   % matrix to save some time
   load Rmatrix_2048_8192.mat
end
%
disp('Finished creating matrix');

hR = @(x) R*x;
hRt = @(x) R'*x;

% noisy observations
sigma = 0.001; 
b = hR(f) + sigma*randn(m,1);

% regularization parameter
 lambda = 0.005*max(abs(R'*b));
%lambdb = 0.05*max(abs(R'*b));

disp('Starting TSIDFPM')
[x1,mses1,Tcpu1,NF1,NormF1] = TSIDFPM(R,b,lambda,f,'TSIDFPM',2,para1);
T1 = Tcpu1(end);

disp('Starting DIDFPM')
[x2,mses2,Tcpu2,NF2,NormF2] =  DIDFPM(R,b,lambda,f,'DIDFPM',2,para2);
T2 = Tcpu2(end);

disp('Starting MTT-CGPM')
[x3,mses3,Tcpu3,NF3,NormF3] = MTTCGPM(R,b,lambda,f,'MTT-CGPM',4,para3);
T3 = Tcpu3(end);

disp('Starting iRMILrp-CGPM')
[x4,mses4,Tcpu4,NF4,NormF4] = iRMILrpCGPM(R,b,lambda,f,'iRMILrp-CGPM',4,para4);
T4 = Tcpu4(end);

disp('Starting PDY')
[x5,mses5,Tcpu5,NF5,NormF5] = PDY(R,b,lambda,f,'PDY',2,para5);
T5 = Tcpu5(end);

fprintf(1,'\n\n-------------------------------------------------\n')   
fprintf(1,'-------------------------------------------------\n')   
fprintf(1,'Problem: n = %g,  m = %g, number of spikes = %g, lambda = %g\n',n,m,n_spikes,lambda)
fprintf(1,'All algorithms initialized with Atb\n')
fprintf(1,'-------------------------------------------------\n')

fprintf(1,'\n TSIDFPMTcpu: %6.2f secs (%d iterations), MSE of the solution = %6.5e\n',...
        T1,length(mses1),mses1(end))  
fprintf(1,'\n DIDFPM Tcpu: %6.2f secs (%d iterations), MSE of the solution = %6.5e\n',...
        T2,length(mses2),mses2(end))
fprintf(1,'\n MTT-CGPM Tcpu: %6.2f secs (%d iterations), MSE of the solution = %6.5e\n',...
        T3,length(mses3),mses3(end))   
fprintf(1,'\n iRMILrp-CGPM Tcpu: %6.2f secs (%d iterations), MSE of the solution = %6.5e\n',...
        T4,length(mses4),mses4(end))
fprintf(1,'\n PDY Tcpu: %6.2f secs (%d iterations), MSE of the solution = %6.5e\n',...
        T5,length(mses5),mses5(end))    
      
fprintf(1,'-------------------------------------------------\n')
fprintf(1,'-------------------------------------------------\n')

% ================= Plotting results =================
%%  semilogy
figure(1)
plot(mses1,'b-','LineWidth',2)
hold on
plot(mses2,'r-','LineWidth',2)
hold on
plot(mses3,'g-','LineWidth',2)
hold on
plot(mses4,'k-','LineWidth',2)
hold on
plot(mses5,'y-','LineWidth',2)
legend('DIDFPM','TSIDFPM','MTT-CGPM','iRMILrp-CGPM','PDY') 
set(gca,'FontName','Times','FontSize',16)
xlabel('Itr')
ylabel('MSE')
title(sprintf('m=%d, n=%d, N=%d, lambda=%.2e',m,n,n_spikes,lambda))
axis on
grid on
hold off

figure(2)
plot(Tcpu1,mses1,'b-','LineWidth',2)
hold on
plot(Tcpu2,mses2,'r-','LineWidth',2)
hold on
plot(Tcpu3,mses3,'g-','LineWidth',2)
hold on
plot(Tcpu4,mses4,'k-','LineWidth',2)
hold on
plot(Tcpu5,mses5,'y-','LineWidth',2)
legend('DIDFPM','TSIDFPM','MTT-CGPM','iRMILrp-CGPM','PDY')
set(gca,'FontName','Times','FontSize',16)
xlabel('Tcpu')
ylabel('MSE')
title(sprintf('m=%d, n=%d, N=%d, lambda=%.2e',m,n,n_spikes,lambda))
axis on
grid on
hold off
%==========================================================================
%%
figure(3)
scrsz = get(0,'ScreenSize');
subplot(6,1,1)
plot(f,'LineWidth',1.1)
top = max(f(:));
bottom = min(f(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original (m = %g, n = %g, N = %g)',m,n,n_spikes))
axis(v)

subplot(6,1,2)
plot(x1(:),'LineWidth',1.1)
top = max(x1(:));
bottom = min(x1(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('DIDFPM(MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',  mses1(end),length(mses1),Tcpu1(end)))
axis(v)

subplot(6,1,3)
plot(x2(:),'LineWidth',1.1)
top = max(x2(:));
bottom = min(x2(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('TSIDFPM(MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',mses2(end),length(mses2),Tcpu2(end)))
axis(v)

subplot(6,1,4)
plot(x3(:),'LineWidth',1.1)
top = max(x3(:));
bottom = min(x3(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('MTT-CGPM(MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',mses3(end),length(mses3),Tcpu3(end)))
axis(v)

subplot(6,1,5)
plot(x4(:),'LineWidth',1.1)
top = max(x4(:));
bottom = min(x4(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('iRMILrp-CGPM (MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',mses4(end),length(mses4),Tcpu4(end)))
axis(v)

subplot(6,1,6)
plot(x5(:),'LineWidth',1.1)
top = max(x5(:));
bottom = min(x5(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('PDY(MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',mses5(end),length(mses5),Tcpu5(end)))
axis(v)
