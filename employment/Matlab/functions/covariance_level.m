function s = covariance_pss_level(t_vector,a_vector,rho_u,rho_v,rho_w,sigma2_e,sigma2_u,sigma2_v,sigma2_theta,sigma2_z)

s=zeros(size(t_vector,1),1);


is_variance=find((a_vector==0));
is_covariance= find((a_vector~=0));


% a!=0 terms
t=t_vector(is_covariance);
a=a_vector(is_covariance);

s(is_covariance) = (1-rho_u.^(t+a+1))/(1-rho_u).*(1-rho_u.^(t+1))/(1-rho_u)*sigma2_theta ... 
    + rho_u.^(2*(t+a+1)-a)*sigma2_u+rho_v.^(2*(t+a+1)-a)*sigma2_v ...
    + sigma2_e*rho_w.^a.*(1-rho_w.^(2*t+2))/(1-rho_w^2);
    
    
% a=0 terms

t=t_vector(is_variance);
a=a_vector(is_variance);

    s(is_variance) =(1-rho_u.^(t+a+1))/(1-rho_u).*(1-rho_u.^(t+1))/(1-rho_u)*sigma2_theta ... 
    + rho_u.^(2*(t+a+1)-a)*sigma2_u+rho_v.^(2*(t+a+1)-a)*sigma2_v ...
    + sigma2_e*rho_w.^a.*(1-rho_w.^(2*t+2))/(1-rho_w^2) ...
    + sigma2_z;

end