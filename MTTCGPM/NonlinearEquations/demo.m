
clc;
clear all
close all

% set random number seed
rng(2016)

ITR_max = 2000;
% setup TXT documentq
fid_tex=fopen('mytext.txt','w');
problem_set = [1:42 64:105 127:189 209:229];
% problem_set = [1:42 64:105 127:189 209:229];

% set parameters
np = length(problem_set); % from problem 1 to problem 15 % the number of the test problems
ns = 5;   % the number of the test algorithms
T = zeros(np,ns);
F = zeros(np,ns);
N = zeros(np,ns);

% parameters for PDY
para1.Itr_max = ITR_max;
para1.gamma = 1;         % the initial guess
para1.sigma = 0.01;      % the coefficient of line search
para1.tau = 0.5;         % the compression ratio
para1.rho = 1;
    
% parameters for DIDFPM 
para2.Itr_max = ITR_max;
para2.gamma = 1;         % the initial guess
para2.thetak = 0.3;
para2.alpha1 = 0.28;
para2.beta = 0.0001;
para2.sigma = 0.001;      % the coefficient of line search
para2.tau = 0.3;         %the compression ratio
para2.rho = 1.8;

% parameters for TSIDFPM
para3.Itr_max = ITR_max;
para3.sigma = 0.001;      % the coefficient of line search
para3.tau = 0.3;         % the compression ratio
para3.alpha = 0.1;
para3.phi = -0.14412;

% parameters for MTT-CGPM
para4.Itr_max = ITR_max;
para4.gamma = 0.45;         % the initial guess
para4.sigma = 0.0001;      % the c oefficient of line search
para4.tau = 0.4;         % the compression ratio
para4.rho = 1.9;
para4.mu = 3;
para4.xi= 1;
para4.eat= 2.07;

% parameters for iRMILrp-CGPM
para5.Itr_max = ITR_max;
para5.gamma = 0.6;         % the initial guess
para5.sigma = 0.01;      % the coefficient of line search
para5.tau = 0.6;         % the compression ratio
para5.rho = 1.1;    

for index=1:np
    Num = problem_set(index);
    %     [name,n] = init(Num);
    [name,x0,x0_old,x0_old1] = init(Num);
    [T1,NFF1,NI1,G1] = PDY(Num,'PDY',2,x0,para1);
    [T2,NFF2,NI2,G2] = DIDFPM(Num,'DIDFPM',2,x0,para2);
    [T3,NFF3,NI3,G3] = TSIDFPM (Num,'TSIDFPM',2,x0,para3);
    [T4,NFF4,NI4,G4] = MTTCGPM (Num,'MTT-CGPM',4,x0,x0_old,x0_old1,para4);
    [T5,NFF5,NI5,G5] = iRMILrpCGPM(Num,'iRMILrp-CGPM',4,x0,x0_old,para5);
        fprintf(fid_tex,'%s & %.3f/%d/%d/%.2e & %.3f/%d/%d/%.2e & %.3f/%d/%d/%.2e&  %.3f/%d/%d/%.2e & %.3f/%d/%d/%.2e\n', ...
        name,T1,NFF1,NI1,G1,T2,NFF2,NI2,G2,T3,NFF3,NI3,G3,T4,NFF4,NI4,G4,T5,NFF5,NI5,G5);
    T(index,:) = [T1,T2,T3,T4,T5];
    F(index,:) = [NFF1,NFF2,NFF3,NFF4,NFF5];
    N(index,:) = [NI1,NI2,NI3,NI4,NI5]; 
end

fclose(fid_tex);

%%画图
clf;   %clf删除当前图形窗口中、
%%句柄未被隐藏(即它们的HandleVisibility属性为on)的图形对象
figure(1);
%subplot(2,2,1);
perf(T,'logplot');
%title('时间性能');
%set(gca,'ylim',[0.3,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('PDY','DIDFPM','TSIDFPM','MTT-CGPM','iRMILrp-CGPM','Location','southeast');
% %subplot(2,2,2);
figure(2);
perf(F,'logplot');
%title('目标函数计算性能');
% set(gca,'ylim',[0.1,1]);
xlabel('\tau','Interpreter','tex');                     %% 给x轴加注
ylabel('\rho(\tau)','Interpreter','tex');               %% 给y轴加注
legend('PDY','DIDFPM','TSIDFPM','MTT-CGPM','iRMILrp-CGPM','Location','southeast')
% %subplot(2,2,3);
figure(3);
perf(N,'logplot');
%title('迭代次数性能');
%set(gca,'ylim',[0.5,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('PDY','DIDFPM','TSIDFPM','MTT-CGPM','iRMILrp-CGPM','Location','southeast');
