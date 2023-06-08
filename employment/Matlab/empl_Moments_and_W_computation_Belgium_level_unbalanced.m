% This code:

% 1) Loads data of firm log size from "new_vat_data_reshaped_level.csv"

% 2) Study missing values in the dataset and does some plots on firm demographics
%    [it can be commented without issues for the subsequent code]

% 3) Computes vector C_level (correlation of levels across ages)

% 4) Computes weighting matrix W_boot_level with bootstrap


%% 1) Loading data
clear
cd 'D:\Kevin\Firm_Growth_Project\_Final\Employment\Data'
addpath('D:\Kevin\Firm_Growth_Project\_Final\Employment\Code\Matlab')
data=readmatrix('D:\Kevin\Firm_Growth_Project\_Final\Employment\Data\clean\reshaped\empl_new_vat_data_reshaped_level_unbalanced.csv');
data=data(:,2:end);
%% 2) Study missing values in the dataset and does some plots on firm demographics

% Nf=size(data,1); % number of firms
% 
% T=100; % number of ages considered
% 
% % We study of number of times a firms appears
% 
% appearing_data=~isnan(data); % logical =1 if value for firm-age is observed
% quarters_appearing=sum(appearing_data,2); % number of times each firm is observed
% 
% % We now control if any pure missing values are present in the dataset
% 
% extreme_dates=zeros(Nf,2); % records date of appearance and disappereance for each firm
% 
% for i=1:Nf
%     extreme_dates(i,1)=min(find(appearing_data(i,:)));
%     extreme_dates(i,2)=max(find(appearing_data(i,:)));
% end
% 
% duration=extreme_dates(:,2)-extreme_dates(:,1)+1;
% 
% if sum((duration~=quarters_appearing))==0
%     display("there are no pure missing values")
%     
% else
%     display("there are pure missing values")
% end
% 
% figure
% histogram(quarters_appearing)
% title("Distribution of firm observed life length")
% xlabel("Number of periods observed")
% ylabel("Quantity of firms")
% 
% 
% figure
% histogram(extreme_dates(:,1))
% title("Distribution of birth date")
% xlabel("Age first time observed")
% ylabel("Quantity of firms")
% 
% figure
% histogram(extreme_dates(:,2))
% title("Distribution of death date")
% xlabel("Age last time observed")
% ylabel("Quantity of firms")
% 
% 

%% 3) Computes vector C_lev (autocovariance of levels across ages)

Nf=size(data,1); % number of firms

T=100; % number of ages considered
a_max=21; %maximum horizon for computation of autocovariance, limit imposed by sample length 1995-2017

% Remark: we could have picked a_max=22, but then N_{t,a=22} would have no values for some ages t (especially old ages)

C_size=(T-a_max)*(a_max+1)+(a_max)*(a_max+1)/2; % number of autocorrelation pairs we can compute
% (T-a_max)*(a_max+1) counts the ages for which we can compute the full (a_max+1) autocov (i.e. for t=1,...,79)
% (a_max)*(a_max+1)/2 counts the ages for which we can compute a decreasing number of autocov (i.e. for t=80,...,100)


t_vector=zeros(C_size,1); % stores the age corresponding to each row of C
a_vector=zeros(C_size,1); % stores the autocov horizon associated to each row of C

for j=1:(T-a_max)
t_vector((j-1)*22+1:j*22)=j*ones(22,1);
a_vector((j-1)*22+1:j*22)=[0:21]';
end

initial_row=1+(T-a_max)*(a_max+1);

for j=(T-a_max+1):T
    max_horizon=T-j;
    t_vector(initial_row:initial_row+max_horizon)=j*ones(max_horizon+1,1);
    a_vector(initial_row:initial_row+max_horizon)=[0:max_horizon]';
    initial_row=initial_row+max_horizon+1;
end

% tot_observed_firms will store total number of firms observed in both t and t+a, that is: N_{t,a}
N_observed_firms=zeros(C_size,1); 

% C is the matrix of autocovariances
% entry (t,t') is autocovariance from age t=1...,T, to age t'=t,...,max(t+a_max,T)
C=zeros(C_size,1);

for j=1:C_size
    t=t_vector(j);
    tprime=t_vector(j)+a_vector(j);
 
        observed_firms_t_tprime=(~isnan(data(:,t))& ~isnan(data(:,tprime))); % selects rows (i.e. firms) of data that has non NaN values (i.e. observed firms)

        firms_t=data(observed_firms_t_tprime,t); % level of observed firms at age t
        firms_tprime=data(observed_firms_t_tprime,tprime); % level of observed firms at age t'
        N_observed_firms(j)=size(firms_t,1); % number of firms observed both at t and t'

        firms_t=firms_t-mean(firms_t);
        firms_tprime=firms_tprime-mean(firms_tprime);
        
        C(j)=(firms_t'*firms_tprime)/N_observed_firms(j); % covariance formula

end

save('matlab\C_level_unbalanced.mat','C')
%save('data\t_vector.mat','t_vector')
%save('data\a_vector.mat','a_vector')

%% 3a) Computing bootstrap weighting matrix

% This section computes bootstrap estimate of W
% Each iteration computes the vector C for half of the sample of firms

% Number of bootstrap iterations:
N_boot=1000;

% Before proceeding, we create a matrix of size Nf x C_size to store for each firm whether it
% is observed both at both horizons corresponding to the row of C

% Another possibility is to use the function "cov" with option "omitrows", that automatically omits NaN
% The issue is that it takes forever because it has to search for NaN in every loop

observed_firms_t_tprime_matrix=sparse(Nf,C_size);

for j=1:C_size
j
   t=t_vector(j);
    tprime=t_vector(j)+a_vector(j);
 
        observed_firms_t_tprime_matrix(:,j)=(~isnan(data(:,t))& ~isnan(data(:,tprime))); % selects rows (i.e. firms) of data that has non NaN values (i.e. observed firms)
end
    
%% Matrix to accomodate bootstrap estimates of C
C_boot=zeros(C_size,N_boot);

% Launching bootstrap iterations

for b=1:N_boot
    b
    %sample half of the firms
selected_firms=(rand(Nf,1)>0.5);
    
    for j=1:C_size
    t=t_vector(j);
    tprime=t_vector(j)+a_vector(j);

    % selects rows (i.e. firms) of data that has non NaN values (i.e. observed firms) and that have been selected by bootstrap
        observed_firms_t_tprime=(observed_firms_t_tprime_matrix(:,j)~=0&selected_firms); 

        firms_t=data(observed_firms_t_tprime,t); % growth of observed firms at age t
        firms_tprime=data(observed_firms_t_tprime,tprime); % growth of observed firms at age t'
        N_observed_firms(j)=size(firms_t,1); % number of firms observed both at t and t'

        firms_t=firms_t-mean(firms_t);
        firms_tprime=firms_tprime-mean(firms_tprime);
        
        C_boot(j,b)=(firms_t'*firms_tprime)/N_observed_firms(j); % covariance formula

    end
    
end

%% building W_boot

% We identify bootstrap estimates that have non NaN values
% Remark: a NaN could emerge if some values of (t,t_prime) do not appear in the bootstrap selection
% In practice this occurs only rarely for very old firms

ok_data_mat=(isnan(C_boot)==0);
N_data=sum(ok_data_mat,2);

W_boot=zeros(C_size,C_size);

for i=1:C_size
    i
    for j=1:C_size
        % (i,j) is the target position of W
        
        
        if (N_data(i)==N_boot & N_data(j)==N_boot) % if there are no NaN values in any of the N_boot estimates
            
            temp1=C_boot(i,:)-mean(C_boot(i,:));
            temp2=C_boot(j,:)-mean(C_boot(j,:));

            W_boot(i,j)=temp1*temp2'/N_boot;
        else % if there are some NaN values in any of the N_boot estimates, we need to take it into account 

            ok_data=(ok_data_mat(i,:)&ok_data_mat(j,:)); % selects rows (i.e. firms) of data that has non NaN values (i.e. observed firms)
            temp1=C_boot(i,ok_data)-mean(C_boot(i,ok_data));
            temp2=C_boot(j,ok_data)-mean(C_boot(j,ok_data));

            W_boot(i,j)=temp1*temp2'/sum(ok_data);

        end
        
    end
end

save('matlab\W_boot_level_unbalanced.mat','W_boot')


