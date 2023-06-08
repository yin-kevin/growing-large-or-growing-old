%% add paths
clear
clc
cd 'F:\Kevin\Firm_Growth_Project\_Final\Employment\Data'

% parameter location
addpath('F:\Kevin\Firm_Growth_Project\_Final\Sales\Estimation\growth_simplified\halt')
addpath('F:\Kevin\Firm_Growth_Project\_Final\Employment\Estimation\growth_simplified\halt')
rng default;

% import parameters
sale_estimates = readtable("HALT_tcap100_iden_growth_simplified_estimates.csv");
empl_estimates = readtable("HALT_tcap100_iden_employment_growth_simplified_estimates.csv"');

%s_rho = sale_estimates{1,2};
s_rho = 0.25;
s_sigma_e2 = sale_estimates{1,3};
s_sigma_z2 = sale_estimates{1,4};
s_psi = sale_estimates{1,5};

e_rho = empl_estimates{1,2};
e_sigma_e2 = empl_estimates{1,3};
e_sigma_z2 = empl_estimates{1,4};
e_psi = empl_estimates{1,5};

% set hyper parameters
alpha = 0.99;
kappa = [-Inf, -1, -0.5, 0, 0.5, 1];
no_of_firms = 10000;


%% sales simulation

s_E = randn(40,no_of_firms).*sqrt(s_sigma_e2);
s_Z = randn(40,no_of_firms).*sqrt(s_sigma_z2);

% compute w_it in matrix form
s_W = cumsum(s_E);

% draw theta_i and initial u_i0 as vectors
s_theta = randn(1,no_of_firms).*sqrt(alpha*s_psi);
s_u0 = randn(1,no_of_firms).*sqrt((1-alpha)*s_psi);

% compute u_it recursively
s_U = s_u0;

for t = 1:39
    % note the index is +1 than the period t
    s_ut = s_U(t,:).*s_rho + s_theta.*(1-s_rho);
    s_U = [s_U;s_ut];
end

% compute log size matrix S and growth rate matrix G
s_S = s_U + s_W + s_Z;



%% sales levels

for i=1:size(kappa,2)
    
    % define levels matrix
    s_L = s_S;
    
    % survival indicator matrix D
    s_D = double(s_U + s_W >= kappa(i));
    s_D(s_D==0) = NaN;
    
    % remove exiting firms from growth rates
    s_L = s_L.*s_D;
    
    % store first row which gets removed in cumsum
    s_L1 = s_L(1,:);
    
    % make exit permanent (cuts another row from the growth rates matrix)
    s_L = [s_L1;diff(cumsum(s_L))];

    % --- compute covariances ---
    s_covL = NaN(size(s_L,1));
    
    for t = 1:size(s_L,1)
       % calculate average growth rate for those alive at t
       g_bar_t = mean(s_L(t,:),"omitnan");

       for a = 0:37
           if t+a <= size(s_L,1)
               % get array of only those firm g-rates that have survived to t+a
               idx = ~isnan(s_L(t+a,:));
               surv = s_L(:,idx);

               % calculate average growth rate for those alive at t+a
               g_bar_ta = mean(surv(t+a,:));

               % compute covariance
               s_covL(t+a,t) = mean((surv(t,:) - g_bar_t).*(surv(t+a,:) - g_bar_ta));
           end
       end
    end

    % --- plots ---
    figure
    plot(s_covL)
    title('Sales-Level Model Simulation')
    subtitle('kappa = ' + string(kappa(i)) + ', alpha = ' + string(alpha))
    xlabel('Firm Age') 
    ylabel('Autocovariance of Sales Level') 
    saveas(gcf,'F:\Kevin\Firm_Growth_Project\_Final\Sales\Simulation\general_sale_level_k' + string(kappa(i)) + '_a' + string(alpha) + '.png')

end


%% sales growth rates


for i=1:size(kappa,2)

    % redefine growth matrix
    s_G = diff(s_S);
    
    % survival indicator matrix D
    s_D = double(s_U + s_W >= kappa(i));
    s_D(s_D==0) = NaN;

    % drop last row
    s_D(end,:) = []; 

    % remove exiting firms from growth rates
    s_G = s_G.*s_D;
   
    % store first row which gets removed in cumsum
    s_G1 = s_G(1,:);
    
    % make exit permanent (cuts another row from the growth rates matrix)
    s_G = [s_G1;diff(cumsum(s_G))];


    % --- compute covariances ---
    s_covG = NaN(size(s_G,1));
    
    for t = 1:size(s_G,1)
       % calculate average growth rate for those alive at t
       g_bar_t = mean(s_G(t,:),"omitnan");

       for a = 0:37
           if t+a <= size(s_G,1)
               % get array of only those firm g-rates that have survived to t+a
               idx = ~isnan(s_G(t+a,:));
               surv = s_G(:,idx);

               % calculate average growth rate for those alive at t+a
               g_bar_ta = mean(surv(t+a,:));

               % compute covariance
               s_covG(t+a,t) = mean((surv(t,:) - g_bar_t).*(surv(t+a,:) - g_bar_ta));
           end
       end
    end

    % --- plots ---
    figure
    plot(s_covG)
    title('Sales-Growth Model Simulation')
    subtitle('kappa = ' + string(kappa(i)) + ', alpha = ' + string(alpha))
    xlabel('Firm Age') 
    ylabel('Autocovariance of Sales Growth') 
    saveas(gcf,'F:\Kevin\Firm_Growth_Project\_Final\Sales\Simulation\general_sale_growth_k' + string(kappa(i)) + '_a' + string(alpha) + '.png')

end


%% employment simulation

e_E = randn(40,no_of_firms).*sqrt(e_sigma_e2);
e_Z = randn(40,no_of_firms).*sqrt(e_sigma_z2);

% compute w_it in matrix form
e_W = cumsum(e_E);

% draw theta_i and initial u_i0 as vectors
e_theta = randn(1,no_of_firms).*sqrt(alpha*e_psi);
e_u0 = randn(1,no_of_firms).*sqrt((1-alpha)*e_psi);

% compute u_it recursively
e_U = e_u0;

for t = 1:39
    % note the index is +1 than the period t
    e_ut = e_U(t,:).*e_rho + e_theta.*(1-e_rho);
    e_U = [e_U;e_ut];
end

% compute log size matrix S and growth rate matrix G
e_S = e_U + e_W + e_Z;


%% employment levels

for i=1:size(kappa,2)
    
    % define levels matrix
    e_L = e_S;
    
    % survival indicator matrix D
    e_D = double(e_U + e_W >= kappa(i));
    e_D(e_D==0) = NaN;
    
    % remove exiting firms from growth rates
    e_L = e_L.*e_D;
    
    % store first row which gets removed in cumsum
    e_L1 = e_L(1,:);
    
    % make exit permanent (cuts another row from the growth rates matrix)
    e_L = [e_L1;diff(cumsum(e_L))];

    % --- compute covariances ---
    e_covL = NaN(size(e_L,1));
    
    for t = 1:size(e_L,1)
       % calculate average growth rate for those alive at t
       g_bar_t = mean(e_L(t,:),"omitnan");

       for a = 0:37
           if t+a <= size(e_L,1)
               % get array of only those firm g-rates that have survived to t+a
               idx = ~isnan(e_L(t+a,:));
               surv = e_L(:,idx);

               % calculate average growth rate for those alive at t+a
               g_bar_ta = mean(surv(t+a,:));

               % compute covariance
               e_covL(t+a,t) = mean((surv(t,:) - g_bar_t).*(surv(t+a,:) - g_bar_ta));
           end
       end
    end

    % --- plots ---
    figure
    plot(e_covL)
    title('Employ-Level Model Simulation')
    subtitle('kappa = ' + string(kappa(i)) + ', alpha = ' + string(alpha))
    xlabel('Firm Age') 
    ylabel('Autocovariance of Employment Level') 
    saveas(gcf,'F:\Kevin\Firm_Growth_Project\_Final\Employment\Simulation\general_empl_level_k' + string(kappa(i)) + '_a' + string(alpha) + '.png')

end


%% employment growth rates

for i=1:size(kappa,2)
    
    % redefine growth matrix
    e_G = diff(e_S);
    
    % survival indicator matrix D
    e_D = double(e_U + e_W >= kappa(i));
    e_D(e_D==0) = NaN;

    % drop last row
    e_D(end,:) = []; 

    % remove exiting firms from growth rates
    e_G = e_G.*e_D;
    
    % store first row which gets removed in cumsum
    e_G1 = e_G(1,:);
    
    % make exit permanent (cuts another row from the growth rates matrix)
    e_G = [e_G1;diff(cumsum(e_G))];


    % --- compute covariances ---
    e_covG = NaN(size(e_G,1));
    
    for t = 1:size(e_G,1)
       % calculate average growth rate for those alive at t
       g_bar_t = mean(e_G(t,:),"omitnan");

       for a = 0:37
           if t+a <= size(e_G,1)
               % get array of only those firm g-rates that have survived to t+a
               idx = ~isnan(e_G(t+a,:));
               surv = e_G(:,idx);

               % calculate average growth rate for those alive at t+a
               g_bar_ta = mean(surv(t+a,:));

               % compute covariance
               e_covG(t+a,t) = mean((surv(t,:) - g_bar_t).*(surv(t+a,:) - g_bar_ta));
           end
       end
    end

    % --- plots ---
    figure
    plot(e_covG)
    title('Employ-Level Model Simulation')
    subtitle('kappa = ' + string(kappa(i)) + ', alpha = ' + string(alpha))
    xlabel('Firm Age') 
    ylabel('Autocovariance of Employment Growth') 
    saveas(gcf,'F:\Kevin\Firm_Growth_Project\_Final\Employment\Simulation\general_empl_growth_k' + string(kappa(i)) + '_a' + string(alpha) + '.png')

end







