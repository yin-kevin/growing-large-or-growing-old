% This code

% 1) Requires to set path for functions

% 2) Requires to set estimation options

% 3) Loads moments and weighting matrix computed in "moments_and_weighting_computation_level_balanced"

% 4) Estimation by minimum distance using fmincon%% load moments and weights

% 5) Run some plots of the results

% 6) Show model fit


% Remark: compared to the file loss_function_minimization_autocov, this code:
% - uses the simplified version of the model, with only N=4 parameters
% - calls the function lossfunction_pss_autocov_level rather than lossfunction_pss_autocov
% - that function in turn calls the function covariance_pss_level rather than covariance_pss

% covariance_pss contained the covariance formulas in equations(1) and (2) of the estimation file (May2021)
% covariance_pss_simplified contains the covariance formulas in equations (5) and (6) of that file
% covariance_pss_level contains the covariance formula for levels in equation (2) of PSS paper
%% 1) Set path for functions
clear
clc
cd 'D:\Kevin\Firm_Growth_Project\_Final\Employment\Data'

% location of code
addpath('D:\Kevin\Firm_Growth_Project\_Final\Employment\Code\Matlab')

% for using functions
addpath('D:\Kevin\Firm_Growth_Project\_Final\Employment\Code\Matlab\functions')

% for storing results
addpath('D:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation')
rng default;

%% 2) Setting estimation options

lower_bound=0;  %set 0 for lower bound=0, 1 for lower bound=-1


W_bootstrap=1; %set 0 for identity weighting, 1 for bootstrap weighting
tcap=19; %max choice is tcap=100;
% t_cap=19 reflects the cap in PSS data

% set subtract_t_vector=1 to compute theoretical moments for t=0,1,2...19
% rather than t=0,1,2,...20
subtract_t_vector=1;

% choose optimization method
optz_method=1;
% optz_method=1 -> fmincon
% optz_method=2 -> fminunc [remark: with this, we cannpt enforce a lower bound - UNC=unconstrained]
% optz_method=3 -> fminsearch
% optz_method=4 -> hybrid: first fminsearch, then fmincon

% needed to cut C as in PSS (tcap=19 and t_plus_a_cap reaching maximum 19) 
PSS_cut=1;
%% 3) Data loading
N = 100; % number of starting points

% Load the datafile, containing the matrix of phi hat values. The row index is "t" and the column
% index is "a"
load('matlab\C_level_balanced.mat'); %already named C once uploaded
load('matlab\t_vector.mat')
load('matlab\a_vector.mat')
load('matlab\W_boot_level_balanced.mat') %already named W_boot once uploaded

% Normalize Bootstrap Matrix 
W_boot=W_boot/(mean(sum(W_boot)));
% this normalization is required for the comparability of Loss values

clear C_level W_boot_level
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

if PSS_cut==0
% impose cap on t
elements_before_tcap=(t_vector<=tcap);
C_cap=C(elements_before_tcap);
t_vector_cap=t_vector(elements_before_tcap);
a_vector_cap=a_vector(elements_before_tcap);
W_choice_cap=W_choice(elements_before_tcap,elements_before_tcap);
end


if PSS_cut==1
% impose cap on and t+a<=19
elements_before_cap=(t_vector+a_vector<=19);
C_cap=C(elements_before_cap);
t_vector_cap=t_vector(elements_before_cap);
a_vector_cap=a_vector(elements_before_cap);
W_choice_cap=W_choice(elements_before_cap,elements_before_cap);    
end

% Normalize Bootstrap Matrix again
W_choice_cap=W_choice_cap/(mean(sum(W_choice_cap)));

L = @(x) lossfunction_level(x,C_cap,t_vector_cap,a_vector_cap,W_choice_cap);

nparams = 8;

x0 = rand(8,1);
L(x0);

% Define the minimizer options and constraints

%for optz_method=1:4 % use this line if you want to run all four optz_method
    
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
a = [0,0,0,0,0,0,0,0]; % lower bound of randomization: Francesco choice
b = [1,1,1,5,5,5,5,5]; % upper bound of randomization


if lower_bound==0  
    lb=zeros(1,8);
else
        lb=ones(1,8)*(-1);
end

ub = [2,2,2,100,100,100,100,100];
x_guess_mat = a + (b-a).*rand(N,nparams);

x_min = zeros(N,nparams);
Loss_value = ones(N,1)*10000;
exitflag = ones(N,1)*(-100);

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
    {'L','rho_u','rho_v','rho_w','sigma2_e','sigma2_u','sigma2_v','sigma2_theta','sigma2_z','exitflag'});

% minimand
x = table2array(allresults_ordered(1,2:9));

% reproduce sound
sound(sin(1:3000));


%% 5) Some plots of the results

varnames=["\rho_u","\rho_v","\rho_w","\sigma^2_e","\sigma^2_u","\sigma^2_v","\sigma^2_\theta","\sigma^2_z"];

figure
for j=1:8
    subplot(3,3,j)
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
writetable(allresults_ordered(1:5,:), 'D:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation\level_balanced\employment_level_balanced_estimates.csv')

%end % use this line if you want to run all four optz_options together
%%
% re-order and take square roots to enhance comparability with PSS
varnames_PSS=["rho_u","rho_v","rho_w"," sigma_t","sigma_u","sigma_v","sigma_e","sigma_z"];

temp=table2array(allresults_ordered(1,[2,3,4,8,6,7,5,9]));
temp=[temp(1:3),temp([4:end]).^0.5];
temp=array2table(temp,"VariableNames",varnames_PSS);
display("This is the estimate that has to be compared with PSS Table 1")
temp

%% 6A) Empirical moments plotting

% First of all let us plot the empirical moments with different cuts

figure
for tt=0:18
    plot([tt:14+tt],C(t_vector==tt),'b') 
    hold on
end
title("Empirical autocovariances")
xlabel("Age (Max age 18)")


figure
for tt=0:40
    plot([tt:14+tt],C(t_vector==tt),'b') 
    hold on
end
title("Empirical autocovariances")
xlabel("Age (Max age 40)")


figure
for tt=0:78
    plot([tt:14+tt],C(t_vector==tt),'b') 
    hold on
end
title("Empirical autocovariances")
xlabel("Age (Max age 78)")


figure
plot(C(a_vector==0),'b') 
title("Empirical Variances")
xlabel("Age")

%% 6B) Show model fit

% Best estimate for the parameters
rho_u_fit=table2array(allresults_ordered(1,2));
rho_v_fit=table2array(allresults_ordered(1,3));
rho_w_fit=table2array(allresults_ordered(1,4));
sigma2_e_fit=table2array(allresults_ordered(1,5));
sigma2_u_fit=table2array(allresults_ordered(1,6));
sigma2_v_fit=table2array(allresults_ordered(1,7));
sigma2_theta_fit=table2array(allresults_ordered(1,8));
sigma2_z_fit=table2array(allresults_ordered(1,9));

% Model implied moments
model_fit=covariance_level(t_vector_cap,a_vector_cap,rho_u_fit,rho_v_fit,rho_w_fit,sigma2_e_fit,sigma2_u_fit,sigma2_v_fit,sigma2_theta_fit,sigma2_z_fit);


if PSS_cut==0
figure
for tt=0:tcap
  plot([tt:14+tt],model_fit(t_vector_cap==tt),'r')
    hold on
    plot([tt:14+tt],C_cap(t_vector_cap==tt),'b') 
    hold on
end
legend("Model","Data")
xlabel("Firm age")
ylabel("Autocovariance")
title("Model fit, Belgian Data, balanced 20Y panel")

else    
figure
for tt=0:19
  plot([tt:19],model_fit(t_vector_cap==tt),'r')
    hold on
    plot([tt:19],C_cap(t_vector_cap==tt),'b') 
    hold on
end
legend("Model","Data")
xlabel("Firm age")
ylabel("Autocovariance")
title("Model fit: Belgian Data, levels, balanced 20Y panel")
saveas(gcf,'D:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation\level_balanced\ModelFit_employment_level_balanced.png')
end

