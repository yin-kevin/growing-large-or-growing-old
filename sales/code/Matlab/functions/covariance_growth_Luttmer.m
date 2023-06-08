function s = covariance_pss(t_vector,a_vector,in_alpha,in_delta,rho,sigma2_q,sigma2_e,lambda)

s=zeros(size(t_vector,1),1);

is_variance=find((a_vector==0));
is_covariance= find((a_vector~=0));


% a!=0 terms
t=t_vector(is_covariance);
a=a_vector(is_covariance);

% define 
alpha = exp(in_alpha) / (1 + exp(in_alpha));
delta = exp(in_delta) / (1 + exp(in_delta));

s(is_covariance) = alpha * delta.^(t+a) .* (1 - alpha * delta.^t) * lambda^2 ...
    + rho.^a .* (rho.^(2*t) * sigma2_q + ((1-rho.^(2*t))/(1-rho^2)) * sigma2_e);
    
% a=0 terms

t=t_vector(is_variance);
a=a_vector(is_variance);

s(is_variance) = alpha * delta.^t .* (1 - alpha * delta.^t) * lambda^2 ...
    + (rho.^(2*t) * sigma2_q + ((1-rho.^(2*t))/(1-rho^2)) * sigma2_e);
  
end