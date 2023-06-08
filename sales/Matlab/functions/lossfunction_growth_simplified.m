function val= lossfunction_growth_simplified(x,C,t_vector,a_vector,W)
%lossfunction_pss_autocov_level

    rho_u = x(1);
    sigma2_e = x(2);
    sigma2_z = x(3);
    psi = x(4);
    
    
    % theoretical moments
    cov_mat = covariance_growth_simplified(t_vector,a_vector,rho_u,sigma2_e,sigma2_z,psi);
        
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