% This code

% 1) Requires to set path for functions

% 2) Requires to set estimation option about lower bound for estimation and weighting matrix

% 3) Loads moments and weighting matrix computed in "moments_and_weighting_computation"

% 4) Estimation by minimum distance using fmincon%% load moments and weights

% 5) Computes GMM standard errors using procedure described in Giuseppe's September 10 2021 email

% 6) Run some plots of the results

% 7) Show model fit

% Remark: compared to the file loss_function_minimization_autocov, this code:
% - uses the simplified version of the model, with only N=4 parameters
% - calls the function lossfunction_pss_autocov_simplified rather than lossfunction_pss_autocov
% - that function in turn calls the function covariance_pss_simplified rather than covariance_pss

% covariance_pss contained the covariance formulas in equations(1) and (2) of the estimation file (May2021)
% covariance_pss_simplified contains the covariance formulas in equations (5) and (6) of that file

%% 1) Set path for functions
clear
clc
cd 'F:\Kevin\Firm_Growth_Project\_Final\Employment\Data'

% location of code
addpath('F:\Kevin\Firm_Growth_Project\_Final\Employment\Code\Matlab')

% for using functions
addpath('F:\Kevin\Firm_Growth_Project\_Final\Employment\Code\Matlab\functions')

% for storing results
addpath('F:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation')
rng default;

%% 2) Setting estimation options

lower_bound=0;  %set 0 for lower bound=0, 1 for lower bound=-1
% Remark: lower bound of zero applies only to sigma_z, all other remain
% with lower bound of -1, but it won't bind

W_bootstrap=0; %set 0 for identity weighting, 1 for bootstrap weighting
tcap=40; %max choice is tcap=100;


% set subtract_t_vector=1 to compute theoretical moments for t=0,1,2...19
% rather than t=0,1,2,...20
subtract_t_vector=0;

% choose optimization method
optz_method=1;
% optz_method=1 -> fmincon
% optz_method=2 -> fminunc [remark: with this, we cannpt enforce a lower bound - UNC=unconstrained]
% optz_method=3 -> fminsearch [remark: with this, we cannpt enforce a lower bound - UNC=unconstrained]
% optz_method=4 -> hybrid: first fminsearch, then fmincon
%% 3) Data loading
N = 200; % number of starting points

% Load the datafile, containing the matrix of phi hat values. The row index is "t" and the column
% index is "a"
load('matlab\C_growth_halt.mat');
load('matlab\t_vector.mat')
load('matlab\a_vector.mat')
load('matlab\W_boot_growth_halt.mat')

% Normalize Bootstrap Matrix
W_boot=W_boot/(mean(sum(W_boot)));

T = 100;
A = 14;

W_identity = eye(size(C,1));

if subtract_t_vector==1
t_vector=t_vector-1;
end
%% 4) Estimation by minimum distance using fmincon
if W_bootstrap==0
    W_choice=W_identity;
else
    W_choice=W_boot;
end

% impose cap on t
elements_before_tcap=(t_vector<=tcap);
C_cap=C(elements_before_tcap);
t_vector_cap=t_vector(elements_before_tcap);
a_vector_cap=a_vector(elements_before_tcap);
W_choice_cap=W_choice(elements_before_tcap,elements_before_tcap);

% Normalize Bootstrap Matrix
W_choice_cap=W_choice_cap/(mean(sum(W_choice_cap)));

L = @(x) lossfunction_growth_simplified(x,C_cap,t_vector_cap,a_vector_cap,W_choice_cap);
%L = @(x) lossfunction_pss_autocov(x,C,T,A,W);
nparams = 4;

x0 = rand(4,1);
L(x0);
% % % Define the minimizer options and constraints
% % fmincon_options = optimoptions(@fmincon,'MaxFunctionEvaluations',50000,...
% %     'StepTolerance',1.0e-6,'FunctionTolerance',1.0e-6,'Algorithm','interior-point','Display','Notify');
% % 


%for optz_method=1:4
    
if optz_method==1
optz_options = optimoptions(@fmincon,'MaxFunctionEvaluations',50000,...
    'StepTolerance',1.0e-6,'FunctionTolerance',1.0e-6,'Algorithm','interior-point','Display','Notify');
end

if optz_method==2 
optz_options = optimoptions(@fminunc,'MaxFunctionEvaluations',50000,...
   'StepTolerance',1.0e-6,'FunctionTolerance',1.0e-6,'Display','Notify');
end

if optz_method==3
optz_options =optimset('TolX',1e-6,'TolFun',1e-6, 'MaxFunEvals', 50000,'MaxIter', 50000,'Display','Notify','LargeScale', 'off');
end

if optz_method==4
optz_options1 =optimset('TolX',1e-6,'TolFun',1e-6, 'MaxFunEvals', 50000,'MaxIter', 50000,'Display','Notify','LargeScale', 'off');

optz_options2 = optimoptions(@fmincon,'MaxFunctionEvaluations',50000,...
    'StepTolerance',1.0e-6,'FunctionTolerance',1.0e-6,'Algorithm','interior-point','Display','Notify');
end



%a = [0,0,0,-1,-1,-1,-1,-1]; % lower bound of randomization: Trevor choice
a = [0,0,0,0]; % lower bound of randomization: Francesco choice
b = [1,5,5,5]; % upper bound of randomization

if lower_bound==0  
    lb = [-1,-1,0,-1]; %make lower bound at zero for variances AND autocorr
else
        lb = [-1,-1,-1,-1];
end

ub = [2,100,100,100];
x_guess_mat = a + (b-a).*rand(N,nparams);

x_min = zeros(N,nparams);
Loss_value = ones(N,1)*10000;
exitflag = ones(N,1)*(-100);

% Define problem
% % problem = createOptimProblem('fmincon','objective',L,'options',fmincon_options,'lb',lb,'ub',ub);
% % 
% % fprintf('Minimizing...\n');
% % tic % Start timing
% % for i = 1:N
% %     problem.x0 = x_guess_mat(i,:);
% %     [x_min(i,:),Loss_value(i),exitflag(i)] = fmincon(problem);
% %     fprintf('%i done\n',i)
% % end



% Define problem
if optz_method==1
problem = createOptimProblem('fmincon','objective',L,'options',optz_options,'lb',lb,'ub',ub);
end

if optz_method==2
problem = createOptimProblem('fminunc','objective',L,'options',optz_options,'lb',lb,'ub',ub);
end

if optz_method==4
problem = createOptimProblem('fmincon','objective',L,'options',optz_options2,'lb',lb,'ub',ub);
end

fprintf('Minimizing...\n');
tic % Start timing
for i = 1:N
    if optz_method==1
            problem.x0 = x_guess_mat(i,:);
    [x_min(i,:),Loss_value(i),exitflag(i)] = fmincon(problem);
    end
    
    if optz_method==2
            problem.x0 = x_guess_mat(i,:);
    [x_min(i,:),Loss_value(i),exitflag(i)] = fminunc(problem);
    end
    
    if optz_method==3
    [x_min(i,:),Loss_value(i),exitflag(i)] = fminsearch(L, x_guess_mat(i,:),optz_options);
    end
    
    if optz_method==4
    [x_min_temp,Loss_value_temp,exitflag_temp] = fminsearch(L, x_guess_mat(i,:),optz_options1);
     problem.x0 = x_min_temp;
    [x_min(i,:),Loss_value(i),exitflag(i)] = fmincon(problem);
    end
    
    fprintf('%i done\n',i);
end





toc
fprintf('Done with minimization\n\n');

allresults = [Loss_value, x_min, exitflag];
allresults_ordered=array2table(sortrows(allresults),'VariableNames',...
    {'L','rho_u','sigma2_e','sigma2_z','Psi','exitflag'});

% minimand
x = table2array(allresults_ordered(1,2:5));

%% 5) GMM Standard error computation
x_top5=table2array(allresults_ordered(1:5,2:5));
std_err_vec=zeros(5,4);
for j=1:5
p=x_top5(j,:);
[cov_mat_p] = GMM_std_errors(t_vector_cap,a_vector_cap,p,W_choice_cap);
std_err_p=diag(cov_mat_p).^(0.5);
std_err_vec(j,:)=std_err_p;
end
%% 6) Some plots of the results

varnames=["rho u","sigma2 e","sigma2 z","Psi"];

figure
for j=1:4
    subplot(2,2,j)
histogram(allresults(:,1+j))
title(varnames(j))
xline(x(j),'--r');
end
sgtitle("Distribution of estimates")


figure
histogram(allresults(:,1))
title("Distribution of Loss")

display("The top five estimates are as follows:")
allresults_ordered(1:5,:)
writetable(allresults_ordered(1:5,:), 'F:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation\growth_simplified\employment_growth_simplified_estimates.csv')

%% 7) Show model fit

% Best estimates from simplified growth model
rho_u_fit=table2array(allresults_ordered(2,2));
sigma2_e_fit=table2array(allresults_ordered(2,3));
sigma2_z_fit=table2array(allresults_ordered(2,4));
psi_fit=table2array(allresults_ordered(2,5));

% Model implied moments
model_fit=covariance_growth_simplified(t_vector_cap,a_vector_cap,rho_u_fit,sigma2_e_fit,sigma2_z_fit,psi_fit);

figure
for tt=1:20
    plot([tt:14+tt],model_fit(t_vector_cap==tt),'r')
    hold on
    plot([tt:14+tt],C_cap(t_vector_cap==tt),'b')    
end
legend("Model","Data")
xlabel("Firm age")
ylabel("Autocovariance")
title("Model fit: Belgian Data, employment growth (simplified model)")
saveas(gcf,'F:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation\growth_simplified\ModelFit_employment_growth_simplified.png')
