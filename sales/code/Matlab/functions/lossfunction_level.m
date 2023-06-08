function val= lossfunction_level(x,C,t_vector,a_vector,W)
%lossfunction_pss_autocov
    rho_u=x(1);
    rho_v=x(2);
    rho_w=x(3);
    sigma2_e=x(4);
    sigma2_u=x(5);
    sigma2_v=x(6);
    sigma2_theta=x(7);
    sigma2_z=x(8);
 
    % theoretical moments
    cov_mat = covariance_level(t_vector,a_vector,rho_u,rho_v,rho_w,sigma2_e,sigma2_u,sigma2_v,sigma2_theta,sigma2_z);
        
    % distance matrix
    e = (cov_mat - C);
    
    %e(isnan(moments(:))) = [];
        val = e'*W*e;
% 
%     if exist('W','var')==1
%         val = e'*W*e;
%     else
%         val = e'*e;
%     end


%     e_phi_w = e_phi.*W_phi;
%     e_var_w = e_var.*W_var;
%     e_w = [e_var_w(:); e_phi_w(:)];
%     e_w(isnan(e_w)) = [];

%    val = sum(e.^2);
end