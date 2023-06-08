function s = covariance_pss(t_vector,a_vector,rho_u,sigma2_e,sigma2_z,psi)

s=zeros(size(t_vector,1),1);


is_variance=find((a_vector==0));
is_covariance= find((a_vector~=0));


% a!=0 terms
t=t_vector(is_covariance);
a=a_vector(is_covariance);

s(is_covariance) = rho_u.^(2*(t-1)) .* rho_u.^a * (1-rho_u)^2 * psi;

s(is_covariance) = s(is_covariance) - sigma2_z.*(a==1);
    
    
% a=0 terms

t=t_vector(is_variance);
a=a_vector(is_variance);

    s(is_variance) = rho_u.^(2*(t-1)) * (1-rho_u)^2 * psi ...
        + sigma2_e + 2*sigma2_z;
  
end